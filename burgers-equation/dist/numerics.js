/* Conservative finite-volume Burgers solvers; usable in browsers and Node. */
(function(root){
'use strict';
const methods={godunov:'Godunov',rusanov:'Rusanov',muscl:'MUSCL · SSP-RK2',maccormack:'MacCormack'};
const cases={
 shock:{name:'Moving shock',left:1,right:0,x0:.3,nu:0,end:.4,bc:'exact'},
 rarefaction:{name:'Rarefaction fan',left:0,right:1,x0:.3,nu:0,end:.4,bc:'exact'},
 entropy:{name:'Transonic rarefaction',left:-1,right:1,x0:.5,nu:0,end:.3,bc:'exact'},
 stationary:{name:'Stationary shock',left:1,right:-1,x0:.5,nu:0,end:.3,bc:'exact'},
 sine:{name:'Smooth wave → shock',nu:0,end:.5,bc:'periodic'},
 viscous:{name:'Viscous traveling shock',left:1,right:0,x0:.3,nu:.01,end:.4,bc:'exact'}
};
const f=u=>.5*u*u;
function problem(key,N=160,nu=.01){return {...cases[key],key,N,dx:1/N,nu:key==='viscous'?nu:0};}
function godunov(a,b){if(a>b)return a+b>=0?f(a):f(b);return a>=0?f(a):b<=0?f(b):0;}
function rusanov(a,b){return .5*(f(a)+f(b))-.5*Math.max(Math.abs(a),Math.abs(b))*(b-a);}
const logcosh=z=>Math.abs(z)+Math.log1p(Math.exp(-2*Math.abs(z)))-Math.LN2;
const shockTime=1/Math.PI;
function referenceAvailable(p,t){return p.key!=='sine'||t<shockTime*(1-1e-10);}
function sinePoint(x,t){
 if(!t)return .5+.5*Math.sin(2*Math.PI*x);
 let lo=x-t,hi=x;
 for(let n=0;n<48;n++){const y=(lo+hi)/2;if(y+t*(.5+.5*Math.sin(2*Math.PI*y))>x)hi=y;else lo=y;}
 return .5+.5*Math.sin(2*Math.PI*(lo+hi)/2);
}
function exactPoint(p,x,t){
 if(p.key==='sine')return referenceAvailable(p,t)?sinePoint(x,t):null;
 const {left:a,right:b,x0}=p,s=.5*(a+b),z=x-x0;
 if(p.nu>0)return s-.5*(a-b)*Math.tanh((a-b)*(z-s*t)/(4*p.nu));
 if(t===0)return z<0?a:z>0?b:.5*(a+b);
 if(a>b)return z<s*t?a:z>s*t?b:.5*(a+b);
 return z<a*t?a:z>b*t?b:z/t;
}
function exactAverage(p,lo,hi,t){
 if(!referenceAvailable(p,t))return null;
 if(p.key==='sine'){
  if(t===0)return .5+.5*(Math.cos(2*Math.PI*lo)-Math.cos(2*Math.PI*hi))/(2*Math.PI*(hi-lo));
  const q=[-.8611363115940526,-.3399810435848563,.3399810435848563,.8611363115940526],w=[.3478548451374538,.6521451548625461,.6521451548625461,.3478548451374538];
  return q.reduce((s,z,i)=>s+.5*w[i]*sinePoint((lo+hi)/2+z*(hi-lo)/2,t),0);
 }
 const {left:a,right:b,x0}=p,s=.5*(a+b),pos=x0+s*t;
 if(p.nu>0){const k=(a-b)/(4*p.nu);return s-.5*(a-b)*(logcosh(k*(hi-pos))-logcosh(k*(lo-pos)))/(k*(hi-lo));}
 if(t===0||a>b)return b+(a-b)*Math.max(0,Math.min(1,((t===0?x0:pos)-lo)/(hi-lo)));
 const P=x=>a*x+(Math.max(x-x0-a*t,0)**2-Math.max(x-x0-b*t,0)**2)/(2*t);
 return (P(hi)-P(lo))/(hi-lo);
}
function initial(p){return Float64Array.from({length:p.N},(_,i)=>exactAverage(p,i*p.dx,(i+1)*p.dx,0));}
function cell(u,p,i,t){return p.bc==='periodic'?u[(i%p.N+p.N)%p.N]:i<0||i>=p.N?exactAverage(p,i*p.dx,(i+1)*p.dx,t):u[i];}
function minmod(a,b,c){return a>0&&b>0&&c>0?Math.min(a,b,c):a<0&&b<0&&c<0?Math.max(a,b,c):0;}
function rhs(u,p,t,key){
 const N=p.N,h=p.dx,flux=new Float64Array(N+1),out=new Float64Array(N);
 const slope=i=>{const l=cell(u,p,i-1,t),v=cell(u,p,i,t),r=cell(u,p,i+1,t);return minmod(2*(v-l),.5*(r-l),2*(r-v));};
 for(let j=0;j<=N;j++){
  const a=cell(u,p,j-1,t),b=cell(u,p,j,t);
  const adv=key==='forward'?f(b):key==='backward'?f(a):key==='rusanov'?rusanov(a,b):key==='muscl'?godunov(a+.5*slope(j-1),b-.5*slope(j)):godunov(a,b);
  flux[j]=adv-p.nu*(b-a)/h;
 }
 for(let i=0;i<N;i++)out[i]=-(flux[i+1]-flux[i])/h;
 return {out,left:flux[0],right:flux[N]};
}
function step(key,u,p,t,dt){
 const first=rhs(u,p,t,key==='maccormack'?'forward':key),v=Float64Array.from(u,(z,i)=>z+dt*first.out[i]);
 if(key!=='muscl'&&key!=='maccormack')return {u:v,boundary:dt*(first.left-first.right),speed:maxSpeed(u,p,t)};
 const second=rhs(v,p,t+dt,key==='maccormack'?'backward':'muscl');
 return {u:Float64Array.from(u,(z,i)=>.5*(z+v[i]+dt*second.out[i])),boundary:.5*dt*(first.left-first.right+second.left-second.right),speed:Math.max(maxSpeed(u,p,t),maxSpeed(v,p,t+dt))};
}
function stats(u,p,t){
 let mass=0,tv=0,min=Infinity,max=-Infinity,l1=0,sq=0,linf=0;const has=referenceAvailable(p,t),errors=new Float64Array(p.N);
 for(let i=0;i<p.N;i++){
  const v=u[i];mass+=v*p.dx;min=Math.min(min,v);max=Math.max(max,v);
  if(i)tv+=Math.abs(v-u[i-1]);
  if(has){const e=v-exactAverage(p,i*p.dx,(i+1)*p.dx,t);errors[i]=e;l1+=Math.abs(e)*p.dx;sq+=e*e*p.dx;linf=Math.max(linf,Math.abs(e));}
 }
 if(p.bc==='periodic')tv+=Math.abs(u[0]-u[p.N-1]);
 return {mass,tv,min,max,l1:has?l1:null,rms:has?Math.sqrt(sq):null,linf:has?linf:null,errors:has?errors:null};
}
function maxSpeed(u,p,t){let m=0;for(const z of u)m=Math.max(m,Math.abs(z));for(const i of [-1,p.N])m=Math.max(m,Math.abs(cell(u,p,i,t)));return m;}
function dtFor(cfl,speed,p){return cfl*p.dx/Math.max(1e-12,speed+2*p.nu/p.dx);}
const API={methods,cases,f,problem,godunov,rusanov,shockTime,referenceAvailable,exactPoint,exactAverage,initial,cell,rhs,step,stats,maxSpeed,dtFor};
if(typeof module!=='undefined'&&module.exports)module.exports=API;else root.Burgers=API;
})(typeof globalThis!=='undefined'?globalThis:this);
