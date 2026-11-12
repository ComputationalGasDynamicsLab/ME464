// Periodic linear convection on [0,2); C=a*dt/dx is signed.
(function(root){
const names={ftbs:'FTBS',ftcs:'FTCS',ftfs:'FTFS',leapfrog:'Leapfrog',btcs:'Implicit BTCS',lw:'Lax–Wendroff',mac:'MacCormack'};
function factors(id,c,t){const s=Math.sin(t),d=1-Math.cos(t);let r=1,i=-c*s;
 if(id==='ftbs')r=1-c*d;
 if(id==='lw'||id==='mac')r=1-c*c*d;
 if(id==='ftfs')r=1+c*d;
 if(id==='btcs'){const z=1+c*c*s*s;return [[1/z,-c*s/z]];}
 if(id==='leapfrog'){const q=c*s,b=1-q*q;if(b>=0)return [[Math.sqrt(b),-q],[-Math.sqrt(b),-q]];return [[0,-q+Math.sqrt(-b)],[0,-q-Math.sqrt(-b)]];}
 return [[r,i]];
}
function status(id,c){if(id==='btcs'||id==='ftfs')return '';if(id==='ftcs')return c===0?'stable':'unstable';if(id==='ftbs')return c>=0&&c<=1?'stable':'unstable';if(id==='leapfrog'&&Math.abs(c)===1)return 'marginal';return Math.abs(c)<=1?'stable':'unstable';}
function step(id,u,prev,c){const N=u.length,v=new Float64Array(N);if(id==='btcs'){
 // Solve the periodic tridiagonal system in O(N) with a rank-one correction.
 const q=c/2,gamma=-(1+q*q),upper=new Float64Array(N),pivot=new Float64Array(N);
 pivot[0]=1-gamma;upper[0]=q/pivot[0];
 for(let j=1;j<N;j++){const diagonal=j===N-1?1-q*q/(1+q*q):1;pivot[j]=diagonal+q*upper[j-1];upper[j]=q/pivot[j];}
 function solve(rhs){const out=new Float64Array(N);out[0]=rhs[0]/pivot[0];for(let j=1;j<N;j++)out[j]=(rhs[j]+q*out[j-1])/pivot[j];for(let j=N-2;j>=0;j--)out[j]-=upper[j]*out[j+1];return out;}
 const rhs=new Float64Array(N);rhs[0]=gamma;rhs[N-1]=q;
 const y=solve(u),z=solve(rhs),corner=-q/gamma,correction=(y[0]+corner*y[N-1])/(1+z[0]+corner*z[N-1]);
 for(let j=0;j<N;j++)v[j]=y[j]-correction*z[j];return v;
 }
 let p;if(id==='mac'){p=new Float64Array(N);for(let j=0;j<N;j++)p[j]=u[j]-c*(u[(j+1)%N]-u[j]);}
 for(let j=0;j<N;j++){const l=(j+N-1)%N,r=(j+1)%N;
 if(id==='ftbs')v[j]=u[j]-c*(u[j]-u[l]);
 if(id==='ftfs')v[j]=u[j]-c*(u[r]-u[j]);
 if(id==='ftcs')v[j]=u[j]-.5*c*(u[r]-u[l]);
 if(id==='lw'||id==='leapfrog'&&!prev)v[j]=u[j]-.5*c*(u[r]-u[l])+.5*c*c*(u[r]-2*u[j]+u[l]);
 if(id==='leapfrog'&&prev)v[j]=prev[j]-c*(u[r]-u[l]);
 if(id==='mac')v[j]=.5*(u[j]+p[j]-c*(p[j]-p[l]));
 }return v;}
root.Convection={names,factors,status,step};
})(typeof module==='object'?module.exports:globalThis);
