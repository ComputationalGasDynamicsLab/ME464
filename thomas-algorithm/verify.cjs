'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm');
const html=fs.readFileSync(__dirname+'/index.html','utf8'),script=[...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m=>m[1]).find(s=>s.includes('const PRESETS'));
new vm.Script(script);
const sandbox={document:{getElementById:()=>({})}};vm.createContext(sandbox);
const kernel=script.slice(script.indexOf('  const $ ='),script.indexOf('  // ---------- render ----------'));
vm.runInContext(kernel+'\nglobalThis.api={load,compute,dominance,factorize,solveFactored,residuals,assembleRod,reuseRodSteps,startRodSession,advanceRodSession,temperatureFigure,checkPrediction,fmt,info,gaussTrace,GAUSS_SYSTEMS,state:()=>S,assign:v=>Object.assign(S,v)};',sandbox);
const A=sandbox.api,close=(a,b,t=2e-11)=>assert.ok(Math.abs(a-b)<t,`${a} vs ${b}`);
function dense(a,b,c,d){
 const n=b.length,M=b.map((v,i)=>Array.from({length:n},(_,j)=>i===j?v:j===i-1?a[i]:j===i+1?c[i]:0)),rhs=Array.from(d);
 for(let j=0;j<n;j++){let q=j;for(let i=j+1;i<n;i++)if(Math.abs(M[i][j])>Math.abs(M[q][j]))q=i;if(M[q][j]===0)throw Error('singular');[M[q],M[j]]=[M[j],M[q]];[rhs[q],rhs[j]]=[rhs[j],rhs[q]];for(let i=j+1;i<n;i++){const f=M[i][j]/M[j][j];for(let k=j;k<n;k++)M[i][k]-=f*M[j][k];rhs[i]-=f*rhs[j];}}
 const x=Array(n).fill(0);for(let i=n-1;i>=0;i--){let v=rhs[i];for(let j=i+1;j<n;j++)v-=M[i][j]*x[j];x[i]=v/M[i][i];}return x;
}
for(const preset of ['warmup','conduction','btcs','convdiff','random','scaled'])for(const n of [3,5,10]){
 A.load(preset,n);const S=A.state(),R=A.compute();assert.equal(R.err,-1);
 const ref=dense(S.a,S.b,S.c,S.d);R.x.forEach((v,i)=>close(v,ref[i]));
 assert.ok(A.residuals(S.a,S.b,S.c,S.d,R.x).scaled<1e-13);assert.equal(A.dominance().ok,true);
 if(preset==='scaled')R.x.forEach(v=>close(v,1));
}
A.load('pivoting',3);assert.equal(A.compute().err,0);dense(A.state().a,A.state().b,A.state().c,A.state().d).forEach(v=>close(v,1));
A.load('singular',3);assert.equal(A.compute().err,2);assert.equal(A.dominance().ok,false);
A.load('tiny',3);let S=A.state(),R=A.compute();assert.equal(R.err,-1);assert.ok(R.relative[0]<1e-10);assert.ok(A.residuals(S.a,S.b,S.c,S.d,R.x).scaled>1e-3);
A.assign({n:3,a:[0,2,1],b:[1,5,2],c:[2,1,0],d:[3,8,3]});assert.equal(A.dominance().ok,true);assert.match(A.dominance().kind,/positive definite/);A.compute().x.forEach(v=>close(v,1));
for(const r of [.01,.5,10]){
 const rod=A.assembleRod(r,100,20,20),F=A.factorize(rod.a,rod.b,rod.c),ref=dense(rod.a,rod.b,rod.c,rod.d);
 A.solveFactored(F,rod.a,rod.d).x.forEach((v,i)=>close(v,ref[i]));
 const demo=A.reuseRodSteps(r,100,20,20);assert.equal(demo.factorizations,1);assert.equal(demo.solves,3);
 let old=[20,20,20];for(const row of demo.rows){const d=old.slice();d[0]+=r*100;d[2]+=r*20;const x=dense(rod.a,rod.b,rod.c,d);row.x.forEach((v,i)=>close(v,x[i]));assert.ok(row.scaled<1e-13);old=Array.from(row.x);}
}
A.load('warmup',3);R=A.compute();assert.equal(A.checkPrediction(R.dp[1].toFixed(4),R,{phase:'fwd',row:1}).correct,true);assert.equal(A.checkPrediction('',R,{phase:'fwd',row:1}).correct,false);assert.equal(A.checkPrediction(-R.dp[1],R,{phase:'fwd',row:1}).correct,false);assert.equal(A.checkPrediction(R.x[1],R,{phase:'back',row:1}).correct,true);
assert.notEqual(A.fmt(1e-14),'0');assert.equal(A.checkPrediction(0,{dp:[1e-14]},{phase:'fwd',row:0}).correct,false);
console.log('PASS: solver versus pivoted reference, scaled systems, pivoting versus singularity, tiny-pivot cancellation, SPD certificate, rod assembly, factor reuse, residuals, and prediction tolerance.');

// Gaussian-elimination prestep: solution, A = LU, no tridiagonal fill-in, Thomas-identical pivots, honest operation counts.
function denseFull(M0,d0){
 const n=d0.length,M=M0.map(r=>r.slice()),rhs=Array.from(d0);
 for(let j=0;j<n;j++){let q=j;for(let i=j+1;i<n;i++)if(Math.abs(M[i][j])>Math.abs(M[q][j]))q=i;if(M[q][j]===0)throw Error('singular');[M[q],M[j]]=[M[j],M[q]];[rhs[q],rhs[j]]=[rhs[j],rhs[q]];for(let i=j+1;i<n;i++){const f=M[i][j]/M[j][j];for(let k=j;k<n;k++)M[i][k]-=f*M[j][k];rhs[i]-=f*rhs[j];}}
 const x=Array(n).fill(0);for(let i=n-1;i>=0;i--){let v=rhs[i];for(let j=i+1;j<n;j++)v-=M[i][j]*x[j];x[i]=v/M[i][i];}return x;
}
for(const key of Object.keys(A.GAUSS_SYSTEMS)){
 const sys=A.GAUSS_SYSTEMS[key],G=A.gaussTrace(key),n=G.n;
 assert.equal(G.err,null);
 denseFull(sys.A,sys.d).forEach((v,i)=>close(v,G.x[i]));
 assert.ok(G.residual<1e-12,`${key} residual ${G.residual}`);
 // A = LU with unit diagonal on L, using the stored multipliers.
 for(let i=0;i<n;i++)for(let j=0;j<n;j++){
  let s=0;for(let k=0;k<n;k++){const l=k===i?1:k<i?(G.L[i][k]||0):0,u=G.U[k][j];if(k<=j)s+=l*u;}
  close(s,sys.A[i][j]);
 }
 // Step list: one step per elimination, one per back-substitution row, plus intro, triangular and summary.
 const elim=G.steps.filter(s=>s.kind==='elim'),solve=G.steps.filter(s=>s.kind==='solve');
 assert.equal(elim.length,n*(n-1)/2);assert.equal(solve.length,n);
 assert.equal(G.steps[0].kind,'intro');assert.equal(G.steps[G.steps.length-1].kind,'summary');
 assert.equal(G.work.elims+G.work.trivial,elim.length);
 assert.ok(G.work.real<=G.work.dense);assert.equal(G.work.dense-G.work.real,G.work.skipped);
 if(sys.kind==='tri'){
  const a=sys.A.map((r,i)=>i?r[i-1]:0),b=sys.A.map((r,i)=>r[i]),c=sys.A.map((r,i)=>i<n-1?r[i+1]:0);
  const F=A.factorize(a,b,c);
  for(let i=0;i<n;i++){
   close(F.m[i],G.U[i][i]);                                   // Thomas pivot is the Gauss pivot
   for(let j=0;j<n;j++){
    if(j>i+1)assert.equal(G.U[i][j],0,`fill-in at U[${i}][${j}] in ${key}`);
    if(j<i-1)assert.equal(G.L[i][j],0,`spurious multiplier L[${i}][${j}] in ${key}`);
   }
  }
  assert.equal(G.work.elims,n-1);assert.equal(G.work.trivial,n*(n-1)/2-(n-1));
  assert.ok(G.work.real<G.work.dense);
  assert.ok(G.work.real<=8*n+2,`${key} used ${G.work.real} operations for n = ${n}`);
 } else {
  assert.equal(G.work.trivial,0);assert.equal(G.work.skipped,0);assert.equal(G.work.real,G.work.dense);
 }
}
A.gaussTrace('dense4').x.forEach((v,i)=>close(v,[1,2,-1,1][i]));
A.gaussTrace('tri4').x.forEach(v=>close(v,1));
A.gaussTrace('tri5').x.forEach(v=>close(v,1));
{const rod=A.assembleRod(.5,100,20,20);A.gaussTrace('rod3').x.forEach((v,i)=>close(v,dense(rod.a,rod.b,rod.c,rod.d)[i]));}
assert.ok(A.gaussTrace('dense4').steps.every(s=>!s.trivial));
assert.ok(A.gaussTrace('tri5').steps.filter(s=>s.trivial).length===6);
console.log('PASS: Gaussian elimination trace reproduces the pivoted solution, satisfies A = LU, creates no tridiagonal fill-in, matches the Thomas pivots, and counts skipped zero work honestly.');

// Full application smoke test with a minimal DOM, including generated prediction fields.
const elements=new Map();function el(id='',attrs=''){const e={id,value:attrs.match(/value="([^"]*)"/)?.[1]||'',checked:false,hidden:false,dataset:{},style:{},listeners:{},classList:{add(){},remove(){},toggle(){}},addEventListener(type,fn){this.listeners[type]=fn;},setAttribute(){},focus(){},scrollIntoView(){}};Object.defineProperty(e,'innerHTML',{get(){return this._html||'';},set(v){this._html=v;parse(v);}});if(id)elements.set(id,e);return e;}
function parse(text){for(const m of text.matchAll(/<[\w-]+\b([^>]*\bid="([^"]+)"[^>]*)>/g))if(!elements.has(m[2]))el(m[2],m[1]);}parse(html);
const rodButtons=[0,1,2].map(i=>Object.assign(el(),{dataset:{assemblyRow:String(i)}}));
const ui={console,document:{getElementById:id=>{assert.ok(elements.has(id),'Missing '+id);return elements.get(id);},querySelectorAll:selector=>selector==='[data-assembly-row]'?rodButtons:[],addEventListener(){},activeElement:null},window:{scrollTo(){}},history:{replaceState(){}},location:{hash:''},setTimeout:fn=>fn(),setInterval:()=>1,clearInterval(){}};
vm.createContext(ui);vm.runInContext(script.replace('  renderAssembly();\n  load("warmup",3);','  globalThis.lab={S,load,compute,render,next,showTab,dGo,dominance,G,gGo,gLoad,gRender,gSteps};\n  renderAssembly();\n  load("warmup",3);'),ui);
assert.equal(ui.lab.S.n,3);elements.get('assembly-load').listeners.click();assert.equal(ui.lab.S.n,3);assert.equal(ui.lab.S.d[0],70);assert.equal(elements.get('der').hidden,true);
elements.get('reuse-run').listeners.click();assert.match(elements.get('reuse-output').innerHTML,/Factorizations: 1/);
ui.lab.load('warmup',3);ui.lab.S.predict=true;ui.lab.S.step=1;ui.lab.S.revealed=false;ui.lab.render();assert.ok(elements.has('prediction-answer'));elements.get('prediction-answer').value=.5;elements.get('step').listeners.click({target:{id:'prediction-check'}});assert.match(elements.get('prediction-feedback').textContent,/Correct/);
ui.lab.S.step=7;ui.lab.S.revealed=true;ui.lab.render();assert.match(elements.get('step').innerHTML,/Scaled residual/);elements.get('check-row').value=2;elements.get('step').listeners.change({target:elements.get('check-row')});
for(const preset of ['pivoting','singular','tiny','scaled']){ui.lab.load(preset);ui.lab.S.step=ui.lab.compute().last;ui.lab.render();}
ui.lab.dGo(10);assert.match(elements.get('d-step').innerHTML,/sufficient|right-hand|remaining|back|pivot/i);
console.log('PASS: real page initializes, assembles/loads the rod, reuses factors, checks predictions, verifies original rows and renders failure examples.');

// The Gauss pane: three tabs switch cleanly, every step of every system renders, and the skip toggle only hides steps.
ui.lab.showTab('gauss');
assert.equal(elements.get('gauss').hidden,false);assert.equal(elements.get('der').hidden,true);assert.equal(elements.get('ex').hidden,true);
assert.match(elements.get('g-step').innerHTML,/Gaussian elimination/);
for(const key of Object.keys(A.GAUSS_SYSTEMS)){
 elements.get('g-system').value=key;elements.get('g-system').listeners.change({target:{value:key}});
 const total=ui.lab.gSteps().length;
 for(let k=0;k<total;k++){
  ui.lab.gGo(k);
  const step=elements.get('g-step').innerHTML,matrix=elements.get('g-mx').innerHTML;
  assert.ok(step.length>80,`${key} step ${k} rendered nothing`);
  assert.doesNotMatch(step,/undefined|NaN/, `${key} step ${k} has a formatting hole`);
  assert.doesNotMatch(matrix,/undefined|NaN/, `${key} matrix ${k} has a formatting hole`);
 }
 assert.match(elements.get('g-step').innerHTML,/Thomas algorithm is/);
 assert.match(elements.get('g-count').innerHTML,/Performed in this run/);
 assert.match(elements.get('g-bridge').innerHTML,/m_i/);
}
elements.get('g-system').listeners.change({target:{value:'tri5'}});
const shown=ui.lab.gSteps().length;
elements.get('g-skip').listeners.change({target:{checked:true}});
const hidden=ui.lab.gSteps().length;
assert.equal(shown-hidden,6,'Hiding zero-multiplier rows should remove exactly the trivial eliminations');
assert.ok(ui.lab.gSteps().every(s=>!s.trivial));
ui.lab.gGo(99);assert.match(elements.get('g-step').innerHTML,/Thomas/);
assert.match(elements.get('g-count').innerHTML,/Spent on known zeros/);
elements.get('g-skip').listeners.change({target:{checked:false}});
elements.get('g-track').listeners.click({target:{dataset:{k:'3'}}});assert.equal(ui.lab.G.step,3);
elements.get('g-next').listeners.click();assert.equal(ui.lab.G.step,4);
elements.get('g-prev').listeners.click();elements.get('g-first').listeners.click();assert.equal(ui.lab.G.step,0);
elements.get('g-last').listeners.click();assert.equal(ui.lab.G.step,ui.lab.gSteps().length-1);
ui.lab.showTab('derive');assert.equal(elements.get('gauss').hidden,true);assert.equal(elements.get('der').hidden,false);
ui.lab.showTab('example');assert.equal(elements.get('gauss').hidden,true);assert.equal(elements.get('ex').hidden,false);
console.log('PASS: Gauss pane renders every step of every system, hides only zero-multiplier eliminations, and the three tabs switch without leaking panes.');

// Notation: the right-hand side is f everywhere the student can see it, never d.
const visible=html.replace(/<script>[\s\S]*?<\/script>/g,'');
for(const stray of ['<var>d</var>','[ A | d ]','d\u2032','d\u2081','\\mathbf{d}'])assert.ok(!visible.includes(stray),`Right-hand side still shown as ${stray}`);
for(const id of ['g-step','g-code','code','step','assembly-reduced','assembly-boundary'])assert.doesNotMatch(String(elements.get(id)._html||elements.get(id).textContent||''),/\bd\u2032|\bd\u2081/,`Right-hand side still shown as d in ${id}`);
assert.match(elements.get('g-code').innerHTML,/f\[i\]/);
console.log('PASS: the right-hand side is written f in every rendered panel, matching the lecture notes.');

// Time-level metadata and figures must use the previous solution, including nonuniform later states.
const session=A.startRodSession({r:.5,left:100,right:20,old:20});
const step1=A.advanceRodSession(session),step2=A.advanceRodSession(session);
close(step1.old[0],20);close(step1.d[0],70);close(step1.time,3.125);
step1.x.forEach((v,i)=>close(step2.old[i],v));close(step2.d[0],step1.x[0]+50);close(step2.time,6.25);
assert.match(A.temperatureFigure(session,step2),/BTCS time step 2/);assert.match(A.temperatureFigure(session,step2),/stroke-dasharray/);
assert.match(elements.get('assembly-old-explanation').textContent,/20.*OLD interior temperature/);
assert.match(elements.get('assembly-equation').textContent,/T1ᵏ = 20/);assert.match(elements.get('assembly-reduced').textContent,/f₁ = 70/);
assert.match(elements.get('reuse-figures').innerHTML,/step 3/);
rodButtons[1].listeners.click();assert.match(elements.get('reuse-figures').innerHTML,/step 3/);
elements.get('reuse-figures').listeners.click({target:{dataset:{rodStep:'1'}}});assert.match(elements.get('reuse-step-detail').innerHTML,/Time step 1/);
elements.get('reuse-next').listeners.click();assert.match(elements.get('reuse-figures').innerHTML,/step 4/);assert.match(elements.get('reuse-step-detail').innerHTML,/Time step 4/);
elements.get('reuse-reset').listeners.click();assert.doesNotMatch(elements.get('reuse-figures').innerHTML,/step 1/);
for(let i=0;i<22;i++)elements.get('reuse-next').listeners.click();assert.equal(elements.get('reuse-next').disabled,true);assert.doesNotMatch(elements.get('reuse-figures').innerHTML,/step 21/);
elements.get('assembly-old').value=32;elements.get('assembly-old').listeners.change();assert.doesNotMatch(elements.get('reuse-figures').innerHTML,/step 1/);assert.match(elements.get('assembly-old-explanation').textContent,/32.*OLD interior temperature/);
console.log('PASS: explicit old-temperature versus boundary RHS, per-time-step figures, correct previous-state reuse, history inspection, restart and parameter-reset behavior.');

// Variable grids retain boundary assembly and fixed physical time; validate a smooth eigenmode.
let previousError;
for(const N of [20,40,80,400]){
 const dt=1e-5,r=.01*dt*N*N,fine=A.startRodSession({N,r,left:0,right:0,old:0});
 fine.history[0].x=Array.from({length:N-1},(_,i)=>Math.sin(Math.PI*(i+1)/N));
 const row=A.advanceRodSession(fine),gain=1/(1+4*r*Math.sin(Math.PI/(2*N))**2);
 assert.equal(row.x.length,N-1);close(row.time,dt);
 row.x.forEach((v,i)=>close(v,gain*Math.sin(Math.PI*(i+1)/N)));
 assert.ok(row.scaled<1e-13);
 const error=Math.abs(row.x[N/2-1]-Math.exp(-.01*Math.PI**2*dt));
 if(previousError&&N<=80)assert.ok(previousError/error>3.8,'Spatial refinement should approach second-order accuracy');
 previousError=error;
}
for(const dt of [3.125,100]){
 const N=400,fine=A.startRodSession({N,r:.01*dt*N*N,left:100,right:20,old:20});
 for(let k=0;k<3;k++){
  const row=A.advanceRodSession(fine);assert.equal(row.x.length,399);
  assert.ok(row.x.every(v=>v>=20-1e-8&&v<=100+1e-8));assert.ok(row.scaled<1e-13);close(row.time,(k+1)*dt);
 }
}
elements.get('profile-n').value=100;elements.get('profile-n').listeners.change();
assert.match(elements.get('reuse-time-note').textContent,/N = 100 intervals, 99 interior unknowns/);
assert.match(elements.get('reuse-time-note').textContent,/Δt = 3.125/);
assert.doesNotMatch(elements.get('reuse-figures').innerHTML,/step 1/);
elements.get('reuse-next').listeners.click();assert.match(elements.get('reuse-step-detail').innerHTML,/99 entries/);
assert.match(elements.get('reuse-time-note').textContent,/r = αΔt\/Δx² = 312.5/);
elements.get('assembly-load').listeners.click();assert.equal(ui.lab.S.n,3);
elements.get('profile-dt').value=.25;elements.get('profile-dt').listeners.change();
assert.match(elements.get('reuse-time-note').textContent,/Δt = 0.25/);assert.doesNotMatch(elements.get('reuse-figures').innerHTML,/step 1/);
elements.get('reuse-run').listeners.click();assert.match(elements.get('reuse-figures').innerHTML,/t = 0.75/);
elements.get('profile-n').value=5000;elements.get('profile-n').listeners.change();assert.equal(elements.get('profile-n').value,400);
console.log('PASS: variable grids through N=400, smooth-mode spatial convergence, large-r maximum principle, independent algebra grid, fixed-time comparisons and input reset/clamping.');
