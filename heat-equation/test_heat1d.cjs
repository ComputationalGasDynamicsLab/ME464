// Test the actual numerical kernels embedded in the teaching page.
// Run with: node test_heat1d.cjs (no dependencies required).
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const html = fs.readFileSync(`${__dirname}/index.html`, 'utf8');
const script = html.match(/<script>([\s\S]*?)<\/script>/)[1];
new vm.Script(script); // Parse the entire application, including the controls.
const kernel = script.slice(script.indexOf('  const L ='), script.indexOf('  const css ='));
const sandbox = { performance, document: { getElementById: () => ({ getContext: () => ({}) }) } };
vm.createContext(sandbox);
vm.runInContext(kernel + `
  globalThis.api = { makeGrid, stepScheme,
    boundary(kind, left = 0, right = 0) { bc = kind; vL = left; vR = right; }
  };`, sandbox);
const { makeGrid, stepScheme, boundary } = sandbox.api;
const close = (a, b, tol = 2e-12) => assert.ok(Math.abs(a-b) < tol, `${a} versus ${b}`);
const sine = (N, mode = 1) => Float64Array.from({length:N+1}, (_,i) =>
  i === 0 || i === N ? 0 : Math.sin(mode*Math.PI*i/N));

// Exact discrete eigenmode amplification checks both matrix assembly and Thomas solve.
boundary('dd');
for (const r of [0.1, 0.4, 0.6, 5]) {
  for (const key of ['ftcs', 'btcs', 'cn']) {
    const N = 10, g = makeGrid(N, r), T = sine(N), q = 4*r*Math.sin(Math.PI/(2*N))**2;
    const amp = key === 'ftcs' ? 1-q : key === 'btcs' ? 1/(1+q) : (1-q/2)/(1+q/2);
    for (let n=0; n<2; n++) assert.ok(stepScheme(key,T,g));
    for (let i=0; i<=N; i++) close(T[i], sine(N)[i]*amp**2);
    close(T[0],0); close(T[N],0);
  }
}
const E = sine(10), I = sine(10);
stepScheme('ftcs', E, makeGrid(10,0.1));
stepScheme('btcs', I, makeGrid(10,0.1));
close(E[5], 0.9902113032590307);
close(I[5], 1/(1+0.4*Math.sin(Math.PI/20)**2));
assert.equal(I[5].toFixed(6),'0.990306');

// The shortest Dirichlet mode grows above the FTCS limit and decays with BTCS.
const unstable = sine(20,19), stable = sine(20,19);
for (let n=0; n<20; n++) {
  stepScheme('ftcs',unstable,makeGrid(20,0.6));
  stepScheme('btcs',stable,makeGrid(20,0.6));
}
assert.ok(Math.max(...unstable.map(Math.abs)) > 100);
assert.ok(Math.max(...stable.map(Math.abs)) < 1e-9);

// Nonzero Dirichlet and mixed-boundary steady solutions must stay steady.
for (const [bc, v0, v1, f] of [
  ['dd',2,2, () => 2], ['dn',1,2, x => 1+2*x], ['nd',2,3,x => 1+2*x]
]) {
  boundary(bc,v0,v1);
  for (const key of ['ftcs','btcs','cn']) {
    const T = Float64Array.from({length:11},(_,i)=>f(i/10));
    stepScheme(key,T,makeGrid(10,0.4));
    T.forEach((v,i)=>close(v,f(i/10)));
  }
}
// Insulated ends conserve discrete total heat (trapezoidal weights).
boundary('nn');
const heat = T => T.reduce((s,v,i)=>s+v*(i===0 || i===T.length-1 ? 0.5 : 1),0);
for (const key of ['ftcs','btcs','cn']) {
  const T = Float64Array.from({length:21},(_,i)=>Math.exp(-80*(i/20-0.35)**2));
  const h = heat(T), g = makeGrid(20,0.4);
  for(let n=0;n<50;n++) stepScheme(key,T,g);
  close(heat(T),h,1e-11);
}
// Fixed-r joint refinement should approach second order for both methods.
boundary('dd');
for (const key of ['ftcs','btcs']) {
  const errors = [10,20,40].map(N => {
    const T = sine(N), g = makeGrid(N,0.25), steps = N*N/25;
    for(let n=0;n<steps;n++) stepScheme(key,T,g);
    return Math.abs(T[N/2] - Math.exp(-0.01*Math.PI**2));
  });
  for(let i=1;i<errors.length;i++) assert.ok(errors[i-1]/errors[i] > 3.8 && errors[i-1]/errors[i] < 4.2);
}
// Exercise the actual independent-step synchronization code, including partial steps.
const compareKernel = script.slice(script.indexOf('  function makeComparison('), script.indexOf('  let comparison ='));
vm.runInContext(compareKernel + `
  globalThis.compareAPI = { makeComparison, advanceComparison,
    configure(end) {
      N = 10; dx = 0.1; tEnd = end; bc = 'dd'; vL = vR = 0; ic = 'sine';
      x = Float64Array.from({length:11}, (_,i) => i*dx);
    }
  };`, sandbox);
const {configure, makeComparison, advanceComparison} = sandbox.compareAPI;
configure(10);
const pair = makeComparison(0.4,5);
// Yield without updating the shared time, then finish the same snapshot.
assert.equal(advanceComparison(pair,0),false);
close(pair.time,0);
assert.equal(advanceComparison(pair,Infinity),true);
close(pair.time,5);
close(pair.runs.ftcs.time,pair.runs.btcs.time);
assert.equal(pair.runs.ftcs.steps,13);
assert.equal(pair.runs.btcs.steps,1);
advanceComparison(pair,Infinity);
close(pair.time,10);
assert.equal(pair.runs.ftcs.steps,26);
assert.equal(pair.runs.btcs.steps,2);
const s = Math.sin(Math.PI/20)**2;
close(pair.runs.ftcs.T[5],((1-4*0.4*s)**12*(1-4*0.2*s))**2);
close(pair.runs.btcs.T[5],(1/(1+20*s))**2);
// Noncommensurate final time and either scheme taking the larger step.
for (const [rF,rB] of [[0.4,5],[0.4,0.1],[0.4,0.4]]) {
  configure(7.3);
  const c = makeComparison(rF,rB);
  while(c.time<c.end) {
    advanceComparison(c,Infinity);
    close(c.runs.ftcs.time,c.time); close(c.runs.btcs.time,c.time);
  }
  close(c.time,7.3);
}
// A long physical run remains synchronized and approaches the expected decay.
configure(100);
const longRun = makeComparison(0.4,5);
while(longRun.time<longRun.end) advanceComparison(longRun,Infinity);
close(longRun.time,100);
close(longRun.runs.ftcs.T[5],((1-4*0.4*s)**12*(1-4*0.2*s))**20);
close(longRun.runs.btcs.T[5],(1/(1+20*s))**20);
// Check the teaching UI's default and the three experiment bindings are present.
assert.match(html,/id="on-btcs" checked/);
for(const id of ['stable','unstable','large']) assert.match(html,new RegExp(`data-lesson="${id}"`));
console.log('PASS: application parses; amplification, boundaries, conservation, convergence, lesson controls, and independent time-step synchronization (including t=100).');
