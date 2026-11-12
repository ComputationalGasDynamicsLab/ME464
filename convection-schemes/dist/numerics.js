// Periodic linear convection on [0,2); C=a*dt/dx is signed.
(function(root){
const names={ftbs:'FTBS',ftcs:'FTCS',leapfrog:'Leapfrog',btcs:'Implicit BTCS',lw:'Lax–Wendroff',mac:'MacCormack'};
function factors(id,c,t){const s=Math.sin(t),d=1-Math.cos(t);let r=1,i=-c*s;
 if(id==='ftbs')r=1-c*d;
 if(id==='lw'||id==='mac')r=1-c*c*d;
 if(id==='btcs')return []; // Analytical amplification is reserved for homework.
 if(id==='leapfrog'){const q=c*s,b=1-q*q;if(b>=0)return [[Math.sqrt(b),-q],[-Math.sqrt(b),-q]];return [[0,-q+Math.sqrt(-b)],[0,-q-Math.sqrt(-b)]];}
 return [[r,i]];
}
function status(id,c){if(id==='btcs')return 'homework';if(id==='ftcs')return c===0?'stable':'unstable';if(id==='ftbs')return c>=0&&c<=1?'stable':'unstable';if(id==='leapfrog'&&Math.abs(c)===1)return 'marginal';return Math.abs(c)<=1?'stable':'unstable';}
function step(id,u,prev,c){const N=u.length,v=new Float64Array(N);if(id==='btcs'){
 // Invert the circulant BTCS operator in Fourier space (exact discrete solve).
 const re=new Float64Array(N),im=new Float64Array(N);
 for(let k=0;k<N;k++){let a=0,b=0;for(let j=0;j<N;j++){const t=2*Math.PI*k*j/N;a+=u[j]*Math.cos(t);b-=u[j]*Math.sin(t);}const q=c*Math.sin(2*Math.PI*k/N),den=1+q*q;re[k]=(a+b*q)/den;im[k]=(b-a*q)/den;}
 for(let j=0;j<N;j++){let s=0;for(let k=0;k<N;k++){const t=2*Math.PI*k*j/N;s+=re[k]*Math.cos(t)-im[k]*Math.sin(t);}v[j]=s/N;}return v;
 }
 let p;if(id==='mac'){p=new Float64Array(N);for(let j=0;j<N;j++)p[j]=u[j]-c*(u[(j+1)%N]-u[j]);}
 for(let j=0;j<N;j++){const l=(j+N-1)%N,r=(j+1)%N;
 if(id==='ftbs')v[j]=u[j]-c*(u[j]-u[l]);
 if(id==='ftcs')v[j]=u[j]-.5*c*(u[r]-u[l]);
 if(id==='lw'||id==='leapfrog'&&!prev)v[j]=u[j]-.5*c*(u[r]-u[l])+.5*c*c*(u[r]-2*u[j]+u[l]);
 if(id==='leapfrog'&&prev)v[j]=prev[j]-c*(u[r]-u[l]);
 if(id==='mac')v[j]=.5*(u[j]+p[j]-c*(p[j]-p[l]));
 }return v;}
root.Convection={names,factors,status,step};
})(typeof module==='object'?module.exports:globalThis);
