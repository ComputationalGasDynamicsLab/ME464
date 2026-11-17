'use strict';
const assert=require('node:assert/strict'),B=require('./dist/numerics.js');
const close=(a,b,tol=2e-12)=>assert.ok(Math.abs(a-b)<tol,`${a} != ${b}`);
function run(key,p,until,cfl=.4){let u=B.initial(p),t=0,budget=0;const mass0=B.stats(u,p,0).mass;let steps=0;while(t<until-1e-14){const dt=Math.min(B.dtFor(cfl,B.maxSpeed(u,p,t),p),until-t),next=B.step(key,u,p,t,dt);u=next.u;budget+=next.boundary;t+=dt;if(++steps>100000)throw Error('Run did not finish');}const stats=B.stats(u,p,t);close(stats.mass-mass0-budget,0,2e-11);return {u,t,stats};}
// Exact scalar Riemann flux, including the transonic entropy case.
close(B.godunov(1,0),.5);close(B.godunov(-2,-1),.5);close(B.godunov(-1,1),0);close(B.godunov(1,-1),.5);close(B.rusanov(2,2),2);
for(const id of Object.keys(B.methods)){
 const p=B.problem('sine',80),u=new Float64Array(80).fill(.7),next=B.step(id,u,p,.1,.001);
 next.u.forEach(v=>close(v,.7));close(next.boundary,0);
 // Periodic nonlinear flow must conserve total mass, including viscous fluxes.
 for(const nu of [0,.01]){const q={...p,nu};run(id,q,.12);}
 // Flux accounting with exact, time-dependent ghost cells.
 run(id,B.problem('viscous',80,.02),.12);
}
// Cell averages resolve jumps, fans and thin viscous transitions analytically.
const shock=B.problem('shock',40);close(B.exactAverage(shock,.49,.51,.4),.5);close(B.exactPoint(shock,.49,.4),1);close(B.exactPoint(shock,.51,.4),0);
const fan=B.problem('entropy',40);close(B.exactPoint(fan,.5,.2),0);close(B.exactAverage(fan,.45,.55,.2),0);
for(const key of ['godunov','rusanov','muscl']){
 const r=run(key,B.problem('entropy',160),.2);
 assert.ok(r.u[79]<0&&r.u[79]>-.1);assert.ok(r.u[80]>0&&r.u[80]<.1);
 assert.ok(r.stats.min>=-1-1e-12&&r.stats.max<=1+1e-12);
 const s=run(key,B.problem('shock',160),.4);
 // Shock center should reach x = 0.5; mass supplies an integral position estimate.
 close(s.stats.mass,.5,5e-5);
}
// Check smooth design order, against cell-average references, at a common time.
for(const key of Object.keys(B.methods)){
 const errors=[80,160,320].map(N=>run(key,B.problem('viscous',N,.04),.1).stats.l1);
 const order=Math.log2(errors[1]/errors[2]);
 assert.ok(order>(key==='godunov'||key==='rusanov'?.8:1.7),`${key} order ${order}: ${errors}`);
 console.log(`${B.methods[key]} viscous-shock observed L1 order: ${order.toFixed(3)}`);
}
const sine=B.problem('sine',80);assert.ok(B.referenceAvailable(sine,.3));assert.equal(B.referenceAvailable(sine,.32),false);assert.equal(B.exactAverage(sine,.1,.2,.4),null);assert.equal(B.stats(B.initial(sine),sine,.4).l1,null);
// Exact average agrees with pointwise quadrature across a narrow viscous layer.
const v=B.problem('viscous',80,.001);let integral=0;const M=10000;for(let j=0;j<M;j++)integral+=B.exactPoint(v,.29+.02*(j+.5)/M,0)/M;close(B.exactAverage(v,.29,.31,0),integral,1e-9);
console.log('PASS: exact Riemann fluxes, entropy fan, shock propagation, periodic conservation, boundary-flux budgets, smooth grid convergence, viscous averages and reference validity.');
// Smooth periodic initial conditions: analytic gradient, shock time and characteristic solution.
for(const key of ['sine','sine2','gaussian']){
 const p=B.problem(key,320),half=B.problem(key,320,.01,p.amplitude/2);
 close(half.tShock,2*p.tShock);
 let minSlope=0;
 for(let j=0;j<20000;j++)minSlope=Math.min(minSlope,B.initialSlope(p,j/20000));
 close(-1/minSlope,p.tShock,2e-8);
 for(const x of [.02,.23,.49,.83]){
  close(B.initialPoint(p,x),B.initialPoint(p,x+1));
  const h=1e-6;close((B.initialPoint(p,x+h)-B.initialPoint(p,x-h))/(2*h),B.initialSlope(p,x),2e-8);
  const t=.7*p.tShock,u=B.exactPoint(p,x,t);close(u,B.initialPoint(p,x-t*u),2e-12);
 }
 assert.equal(B.referenceAvailable(p,p.tShock),false);
 assert.equal(B.exactAverage(p,.2,.3,1.5*p.tShock),null);
 const u0=B.initial(p),compress=u=>Math.max(...Array.from(u,(v,i)=>(v-u[(i+1)%p.N])/p.dx));
 const r=run('godunov',p,1.5*p.tShock);
 assert.ok(compress(r.u)>2*compress(u0),`${key}: descending front did not steepen`);
 assert.ok(r.stats.min>=p.range[0]-1e-12&&r.stats.max<=p.range[1]+1e-12);
}
console.log('PASS: smooth profile periodicity, gradients, inverse-amplitude shock times, characteristics, conservation and nonlinear steepening.');
