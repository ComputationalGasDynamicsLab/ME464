'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm');
const html=fs.readFileSync(__dirname+'/index.html','utf8'),script=[...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m=>m[1]).find(s=>s.includes('const PRESETS'));
new vm.Script(script);
const sandbox={document:{getElementById:()=>({})}};vm.createContext(sandbox);
const kernel=script.slice(script.indexOf('  const $ ='),script.indexOf('  // ---------- render ----------'));
vm.runInContext(kernel+'\nglobalThis.api={load,compute,dominance,factorize,solveFactored,residuals,assembleRod,reuseRodSteps,checkPrediction,fmt,info,state:()=>S,assign:v=>Object.assign(S,v)};',sandbox);
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

// Full application smoke test with a minimal DOM, including generated prediction fields.
const elements=new Map();function el(id='',attrs=''){const e={id,value:attrs.match(/value="([^"]*)"/)?.[1]||'',checked:false,hidden:false,dataset:{},style:{},listeners:{},classList:{add(){},remove(){},toggle(){}},addEventListener(type,fn){this.listeners[type]=fn;},setAttribute(){},focus(){},scrollIntoView(){}};Object.defineProperty(e,'innerHTML',{get(){return this._html||'';},set(v){this._html=v;parse(v);}});if(id)elements.set(id,e);return e;}
function parse(text){for(const m of text.matchAll(/<[\w-]+\b([^>]*\bid="([^"]+)"[^>]*)>/g))if(!elements.has(m[2]))el(m[2],m[1]);}parse(html);
const rodButtons=[0,1,2].map(i=>Object.assign(el(),{dataset:{assemblyRow:String(i)}}));
const ui={console,document:{getElementById:id=>{assert.ok(elements.has(id),'Missing '+id);return elements.get(id);},querySelectorAll:selector=>selector==='[data-assembly-row]'?rodButtons:[],addEventListener(){},activeElement:null},window:{scrollTo(){}},history:{replaceState(){}},location:{hash:''},setTimeout:fn=>fn(),setInterval:()=>1,clearInterval(){}};
vm.createContext(ui);vm.runInContext(script.replace('  renderAssembly();\n  load("warmup",3);','  globalThis.lab={S,load,compute,render,next,showTab,dGo,dominance};\n  renderAssembly();\n  load("warmup",3);'),ui);
assert.equal(ui.lab.S.n,3);elements.get('assembly-load').listeners.click();assert.equal(ui.lab.S.n,3);assert.equal(ui.lab.S.d[0],70);assert.equal(elements.get('der').hidden,true);
elements.get('reuse-run').listeners.click();assert.match(elements.get('reuse-output').innerHTML,/Factorizations: 1/);
ui.lab.load('warmup',3);ui.lab.S.predict=true;ui.lab.S.step=1;ui.lab.S.revealed=false;ui.lab.render();assert.ok(elements.has('prediction-answer'));elements.get('prediction-answer').value=.5;elements.get('step').listeners.click({target:{id:'prediction-check'}});assert.match(elements.get('prediction-feedback').textContent,/Correct/);
ui.lab.S.step=7;ui.lab.S.revealed=true;ui.lab.render();assert.match(elements.get('step').innerHTML,/Scaled residual/);elements.get('check-row').value=2;elements.get('step').listeners.change({target:elements.get('check-row')});
for(const preset of ['pivoting','singular','tiny','scaled']){ui.lab.load(preset);ui.lab.S.step=ui.lab.compute().last;ui.lab.render();}
ui.lab.dGo(10);assert.match(elements.get('d-step').innerHTML,/sufficient|right-hand|remaining|back|pivot/i);
console.log('PASS: real page initializes, assembles/loads the rod, reuses factors, checks predictions, verifies original rows and renders failure examples.');
