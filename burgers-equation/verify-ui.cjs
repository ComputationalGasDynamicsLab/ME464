// Run the real app with a minimal DOM/canvas adapter; no browser dependencies.
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),B=require('./dist/numerics.js');
const html=fs.readFileSync(__dirname+'/dist/index.html','utf8'),els=new Map();
const ctx=new Proxy({},{get:(o,k)=>o[k]??(()=>{}),set:(o,k,v)=>(o[k]=v,true)});
function element(id='',attrs=''){
 const e={id,value:attrs.match(/\bvalue="([^"]*)"/)?.[1]||'',checked:/\bchecked\b/.test(attrs),disabled:false,clientWidth:440,clientHeight:id.startsWith('error-')?105:200,style:{setProperty(){}},classList:{add(){},remove(){},toggle(){}},getContext:()=>ctx,appendChild(child){parse(child.innerHTML||'');},click(){},textContent:''};
 if(id)els.set(id,e);return e;
}
function parse(text){for(const m of text.matchAll(/<([a-z]+)\b([^>]*\bid="([^"]+)"[^>]*)>/g)){if(els.has(m[3]))continue;element(m[3],m[2]);}}
parse(html);
for(const m of html.matchAll(/<select[^>]*id="([^"]+)"[^>]*>([\s\S]*?)<\/select>/g)){
 const opts=[...m[2].matchAll(/<option([^>]*)>([^<]*)<\/option>/g)],o=opts.find(z=>/selected/.test(z[1]))||opts[0];els.get(m[1]).value=o[1].match(/value="([^"]*)"/)?.[1]||o[2];
}
const presets=[...html.matchAll(/data-preset="([^"]+)"/g)].map(m=>Object.assign(element(),{dataset:{preset:m[1]}}));
const sandbox={Burgers:B,console,performance,Float64Array,Blob,window:{devicePixelRatio:1},document:{getElementById:id=>{assert.ok(els.has(id),'Missing element '+id);return els.get(id);},createElement:()=>element(),querySelectorAll:()=>presets},ResizeObserver:class{observe(){}},requestAnimationFrame(){},setTimeout(){},URL};
vm.createContext(sandbox);vm.runInContext(fs.readFileSync(__dirname+'/dist/app.js','utf8')+'\nglobalThis.lab={reset,advance,draw,csv,state:()=>({p,t,n,states})};',sandbox);
const lab=sandbox.lab;
assert.equal(lab.state().n,0);els.get('step').onclick();assert.equal(lab.state().n,1);
for(const s of Object.values(lab.state().states))assert.equal(s.time,lab.state().t);
assert.equal(lab.csv().split('\n').length,1+160*4);
// A complete run lands on the final time, with finite errors and matching boundary budgets.
while(lab.advance()){}lab.draw();assert.equal(lab.state().t,.4);
for(const s of Object.values(lab.state().states)){assert.ok(Number.isFinite(s.stats.l1));assert.ok(Math.abs(s.stats.mass-s.initialMass-s.budget)<2e-12);}
// Exact-reference cutoff and CSV must not fabricate post-shock sine errors.
els.get('case').value='sine';els.get('end').value=.34;lab.reset();while(lab.advance()){}lab.draw();for(const s of Object.values(lab.state().states))assert.equal(s.stats.l1,null);
assert.equal(lab.csv().split('\n')[1].split(',')[6],'');
// Unsafe experiment freezes each method independently and keeps its recorded time.
presets.find(z=>z.dataset.preset==='unsafe').onclick();for(let i=0;i<60;i++)lab.advance();lab.draw();
const stopped=Object.values(lab.state().states).filter(s=>s.frozen);assert.ok(stopped.length>0);for(const s of stopped)assert.ok(s.time<=lab.state().t);
// Fit axes, hidden exact reference, viscous controls, and an empty history all render.
els.get('scale').value='fit';lab.draw();presets.find(z=>z.dataset.preset==='viscous').onclick();els.get('step').onclick();lab.draw();assert.equal(lab.state().p.nu,.01);
console.log('PASS: real app initializes, steps, finishes at exact time, exports CSV, suppresses post-shock exact errors, freezes unstable runs, and renders viscous cases.');
