const assert=require('node:assert/strict');
const {Convection:C}=require('./dist/numerics.js');
const N=40,k=7,theta=2*Math.PI*k/N,u=Float64Array.from({length:N},(_,j)=>Math.cos(theta*j));
let checks=0;
function near(a,b,tol=2e-12){assert.ok(Math.abs(a-b)<tol,`${a} vs ${b}`);checks++;}
for(const id of Object.keys(C.names))for(const c of [-1.2,-.8,0,.8,1.2,2]){
 if(id==='leapfrog')continue;
 const [re,im]=C.factors(id,c,theta)[0],v=C.step(id,u,null,c);
 for(let j=0;j<N;j++)near(v[j],re*Math.cos(theta*j)-im*Math.sin(theta*j));
 // Conservative periodic updates preserve the mean.
 near(v.reduce((a,b)=>a+b,0),u.reduce((a,b)=>a+b,0));
}
for(const c of [-1.2,-.8,0,.8,1.2]){const a=C.step('lw',u,null,c),b=C.step('mac',u,null,c);a.forEach((v,j)=>near(v,b[j]));}
// BTCS must satisfy its coupled discrete equation, including the wraparound rows.
for(const c of [-2,.8,2]){const v=C.step('btcs',u,null,c);for(let j=0;j<N;j++)near(v[j]+c/2*(v[(j+1)%N]-v[(j+N-1)%N]),u[j]);}
for(const size of [20,21,81,320,2000])for(const c of [-10,-2,0,.8,2,10]){
 const input=Float64Array.from({length:size},(_,j)=>1+Math.sin(.7*j)+.1*Math.cos(.19*j)),out=C.step('btcs',input,null,c);
 for(let j=0;j<size;j++)near(out[j]+c/2*(out[(j+1)%size]-out[(j+size-1)%size]),input[j],2e-11);
 near(out.reduce((a,b)=>a+b,0)/size,input.reduce((a,b)=>a+b,0)/size,2e-11);
}
// FTFS is exercised through its discrete update, with no stability answer supplied.
for(const c of [-2,-.8,0,.8,2]){const out=C.step('ftfs',u,null,c);for(let j=0;j<N;j++)near(out[j],u[j]-c*(u[(j+1)%N]-u[j]));}
// Leapfrog recurrence, with a deliberately independent older time level.
const old=Float64Array.from(u,v=>.7*v),v=C.step('leapfrog',u,old,.8);for(let j=0;j<N;j++)near(v[j],old[j]-.8*(u[(j+1)%N]-u[(j+N-1)%N]));
for(const [id,c,stable]of [['ftbs',.8,true],['ftbs',1.2,false],['ftcs',.1,false],['lw',1,true],['lw',1.2,false],['mac',-.8,true],['leapfrog',.8,true],['leapfrog',1.2,false]]){let max=0;for(let j=0;j<=1000;j++)for(const z of C.factors(id,c,Math.PI*j/1000))max=Math.max(max,Math.hypot(...z));assert.equal(max<=1+1e-12,stable);checks++;}
// At the leapfrog endpoint a repeated root admits linearly growing two-level data.
let prev=Float64Array.from({length:40},(_,j)=>Math.cos(Math.PI*j/2)),cur=new Float64Array(40);for(let n=0;n<20;n++){let next=C.step('leapfrog',cur,prev,1);prev=cur;cur=next;}assert.ok(Math.max(...cur.map(Math.abs))>10);checks++;
assert.ok(Object.hasOwn(C.names,'ftfs'));assert.equal(C.status('ftfs',.8),'');assert.equal(C.factors('ftfs',.8,theta).length,1);checks+=2;assert.equal(C.status('btcs',2),'');assert.equal(C.factors('btcs',2,theta).length,1);checks+=3;
console.log(`${checks} checks passed: Fourier symbols, BTCS residual, conservation, LW/MacCormack equivalence, stability limits, leapfrog endpoint.`);
