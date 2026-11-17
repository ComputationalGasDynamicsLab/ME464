'use strict';
const B=Burgers,ids=Object.keys(B.methods),$=id=>document.getElementById(id);
const colors={godunov:'#16725c',rusanov:'#4269b2',muscl:'#bf603b',maccormack:'#93657d'};
let p,t=0,n=0,lastDt=0,states={},running=false,lastFrame=0,credit=0,fixedDt=0,end=.4,history=[];
const fmt=x=>x===null?'—':!Number.isFinite(x)?'nonfinite':x===0?'0':Math.abs(x)<.001||Math.abs(x)>=1e4?x.toExponential(2):String(+x.toPrecision(5));
for(const id of ids){const card=document.createElement('article');card.className='panel method-card';card.style.setProperty('--method',colors[id]);card.innerHTML=`<div class="method-head"><h3>${B.methods[id]}</h3><span id="badge-${id}" class="badge"></span></div><p id="time-${id}" class="method-time"></p><div id="metrics-${id}" class="metrics"></div><p class="chart-title">CELL AVERAGES · u(x)</p><canvas id="profile-${id}" role="img" aria-label="${B.methods[id]}: numerical Burgers solution and exact reference"></canvas><p class="chart-title">SIGNED ERROR · independently fitted axes</p><canvas class="error-chart" id="error-${id}" role="img" aria-label="${B.methods[id]}: cell-average error against the exact solution"></canvas><p id="note-${id}" class="card-note"></p>`;$('plots').appendChild(card);}
function updatePrimerProblem(){
 let formula,note;
 if(p.smooth){
  formula=p.key==='gaussian'?`u₀(x) = 0.2 + ${fmt(p.amplitude)} exp{4[cos(2π(x−0.3))−1]}`:`u₀(x) = 0.5 + ${fmt(p.amplitude)} sin(${2*p.waves}πx)`;
  note=`Smooth periodic initial data. For this amplitude, tₛ = ${fmt(p.tShock)}. The falling side compresses; the numerical solutions continue after the first shock forms.`;
 }else if(p.nu>0){formula=`u₀(x) = 0.5 − 0.5 tanh[(x−0.3)/(4ν)], ν = ${fmt(p.nu)}`;note='This initial condition is continuous. It is an exact traveling-layer profile, rather than an inviscid smooth-wave breaking experiment.';}
 else{formula=`u₀(x) = ${p.left} for x < ${p.x0}; ${p.right} for x > ${p.x0}`;note='This Riemann initial condition has a jump. The code integrates the initial profile over each cell, including any cell cut by the jump.';}
 $('primer-ic-formula').textContent=formula;$('primer-ic-note').textContent=note;
 $('primer-periodic').hidden=p.bc!=='periodic';$('primer-exact-boundary').hidden=p.bc==='periodic';
 $('primer-bc-note').textContent=p.bc==='periodic'?'Selected problem: periodic copies at both ends. No independent endpoint value is prescribed.':'Selected problem: exact exterior cell averages at each stage time. The domain can exchange the conserved quantity with its exterior.';
}
function renderPrimerExample(){
 const id=$('primer-scheme').value,example=B.problem('sine',5),u=Float64Array.from([.2,.6,1,.6,.2]),next=B.step(id,u,example,0,.05).u;
 const cells=a=>Array.from(a,v=>`<span>${v.toFixed(4)}</span>`).join('');
 $('primer-old').innerHTML=cells(u);$('primer-new').innerHTML=cells(next);
 const m=a=>Array.from(a).reduce((sum,v)=>sum+v*.2,0);
 $('primer-balance').textContent=`${B.methods[id]}: center ${u[2].toFixed(4)} → ${next[2].toFixed(4)}. Total ΣūᵢΔx: ${m(u).toFixed(6)} → ${m(next).toFixed(6)}. Periodic transfers conserve the total.`;
}

function number(id,lo,hi,fallback){let z=Number($(id).value);if(!Number.isFinite(z))z=fallback;z=Math.min(hi,Math.max(lo,z));$(id).value=z;return z;}
function pause(){running=false;$('play').textContent=t>=end?'Replay':'▶ Play';}
function reset(){
 pause();p=B.problem($('case').value,number('grid',40,400,160),number('viscosity',.001,.1,.01),number('amplitude',.1,1,.5));end=number('end',.01,5,p.end);
 number('cfl',.1,1.5,.45);t=0;n=0;lastDt=0;credit=0;states={};history=[];
 for(const id of ids){const u=B.initial(p),s=B.stats(u,p,0);states[id]={u,time:0,steps:0,initialMass:s.mass,budget:0,stats:s,load:0,frozen:false};}
 fixedDt=B.dtFor(+$('cfl').value,B.maxSpeed(states.godunov.u,p,0),p);
 $('grid-value').textContent=p.N;$('nu-value').textContent=p.nu;$('cfl-value').textContent=fmt(+$('cfl').value);$('viscosity').disabled=p.key!=='viscous';
 $('amplitude').disabled=!p.smooth;$('amplitude-value').textContent=fmt(p.amplitude??+$('amplitude').value);
 $('problem-note').textContent=p.smooth?`${p.name}. Periodic domain [0,1). Predicted shock formation tₛ = ${fmt(p.tShock)}; the simulation continues beyond it.`:`${p.name}: uL = ${p.left}, uR = ${p.right}, x₀ = ${p.x0}. Ghost-cell averages follow the exact whole-line solution at each stage.`;
 updatePrimerProblem();
 $('run-message').textContent='Paused at the initial condition. Predict the wave, then press Play.';$('run-message').className='run-message';draw();
}
function advance(){
 if(t>=end-1e-14||n>=20000){pause();return false;}
 let speed=0;for(const id of ids)if(!states[id].frozen)speed=Math.max(speed,B.maxSpeed(states[id].u,p,t));
 if(ids.every(id=>states[id].frozen)){pause();return false;}
 const requested=$('timestep').value==='fixed'?fixedDt:B.dtFor(+$('cfl').value,speed,p),dt=Math.min(requested,end-t);
 if(!Number.isFinite(dt)||dt<1e-12){pause();$('run-message').textContent='Time step became too small. Reset with a lower factor or a coarser grid.';return false;}
 for(const id of ids){
  const s=states[id];if(s.frozen)continue;
  const next=B.step(id,s.u,p,t,dt);s.u=next.u;s.budget+=next.boundary;s.time=t+dt;s.steps++;s.load=next.speed*dt/p.dx+2*p.nu*dt/p.dx**2;s.stats=B.stats(s.u,p,s.time);
  s.frozen=s.u.some(z=>!Number.isFinite(z)||Math.abs(z)>1e3);
 }
 t+=dt;n++;lastDt=dt;
 const row={time:t};for(const id of ids)row[id]=states[id].frozen?null:states[id].stats.l1;history.push(row);
 if(history.length>1000)history=history.filter((_,i)=>i%2===0);
 if(t>=end-1e-14){t=end;pause();}return true;
}
function chart(id,xrange,yrange,label){
 const cv=$(id),c=cv.getContext('2d'),w=cv.clientWidth,h=cv.clientHeight,dpr=window.devicePixelRatio||1;if(!w||!h)return null;
 cv.width=Math.round(w*dpr);cv.height=Math.round(h*dpr);c.setTransform(dpr,0,0,dpr,0,0);
 const box={l:59,r:w-12,t:20,b:h-28},X=x=>box.l+(x-xrange[0])/(xrange[1]-xrange[0])*(box.r-box.l),Y=y=>box.b-(y-yrange[0])/(yrange[1]-yrange[0])*(box.b-box.t);
 c.font='9px ui-monospace,monospace';c.lineWidth=1;
 for(const a of [0,.5,1]){const x=xrange[0]+a*(xrange[1]-xrange[0]),y=yrange[0]+a*(yrange[1]-yrange[0]);c.strokeStyle='#e0e4d9';c.beginPath();c.moveTo(X(x),box.t);c.lineTo(X(x),box.b);c.stroke();c.beginPath();c.moveTo(box.l,Y(y));c.lineTo(box.r,Y(y));c.stroke();c.fillStyle='#6e796e';c.textAlign='right';c.fillText(fmt(y),box.l-6,Y(y)+3);c.textAlign='center';c.fillText(fmt(x),X(x),h-10);}
 c.strokeStyle='#cbd4c5';c.strokeRect(box.l,box.t,box.r-box.l,box.b-box.t);c.fillStyle='#6e796e';c.textAlign='left';c.fillText(label,box.l,10);
 function line(points,color,dash=[],dots=false){c.save();c.beginPath();c.rect(box.l,box.t,box.r-box.l,box.b-box.t);c.clip();c.strokeStyle=c.fillStyle=color;c.lineWidth=1.7;c.setLineDash(dash);c.beginPath();let started=false;for(const [x,y]of points){if(y===null||!Number.isFinite(y)){started=false;continue;}const yy=Math.max(yrange[0]-20*(yrange[1]-yrange[0]),Math.min(yrange[1]+20*(yrange[1]-yrange[0]),y));if(started)c.lineTo(X(x),Y(yy));else c.moveTo(X(x),Y(yy));started=true;}c.stroke();if(dots&&points.length<=160)for(const [x,y]of points)if(Number.isFinite(y)&&y>=yrange[0]&&y<=yrange[1]){c.beginPath();c.arc(X(x),Y(y),1.7,0,2*Math.PI);c.fill();}c.restore();}
 return {line,c,X,Y,box};
}
function draw(){
 if(!p)return;$('time').textContent=t.toFixed(4);$('step-count').textContent=n;$('dx').textContent=fmt(p.dx);$('dt').textContent=lastDt?fmt(lastDt):'—';
 const init=Array.from({length:401},(_,i)=>[i/400,B.exactPoint(p,i/400,0)]),showInitial=$('initial').checked,cache=new Map();
 const baselo=p.smooth?p.range[0]-.1:Math.min(p.left,p.right)-.15,basehi=p.smooth?p.range[1]+.1:Math.max(p.left,p.right)+.15;
 $('smooth-summary').hidden=!p.smooth;
 if(p.smooth){
  $('initial-formula').textContent=p.key==='gaussian'?`u₀(x) = 0.2 + ${fmt(p.amplitude)} exp{4[cos(2π(x−0.3))−1]}`:`u₀(x) = 0.5 + ${fmt(p.amplitude)} sin(${p.waves===1?'2':'4'}πx)`;
  $('shock-time').textContent=`tₛ = ${fmt(p.tShock)}`;
  $('shock-stage').textContent=t===0?'Smooth initial profile':t<p.tShock?`Steepening · t/tₛ = ${fmt(t/p.tShock)}`:`After breaking · t/tₛ = ${fmt(t/p.tShock)}`;
  const ic=chart('initial-profile',[0,1],[baselo,basehi],'u₀(x) · initial condition');if(ic)ic.line(init,'#16725c');
 }
 const primerPlot=chart('primer-initial',[0,1],[baselo,basehi],'u₀(x) · selected initial condition');if(primerPlot)primerPlot.line(init,'#16725c');
 for(const id of ids){
  const s=states[id],d=s.stats,guide=id==='muscl'?.5:1,unresolved=p.nu>0&&4*p.nu/(p.left-p.right)<3*p.dx;
  $('badge-'+id).textContent=s.frozen?'Frozen':!n?'Ready':s.load>guide+1e-10?'CFL guideline exceeded':id==='maccormack'?'Not monotone':'Within CFL guideline';
  $('badge-'+id).className='badge '+(s.frozen?'bad':s.load>guide+1e-10?'warn':'');
  $('time-'+id).textContent=`t = ${fmt(s.time)} · ${s.steps} steps${s.frozen?' · FROZEN':''}`;
  const residual=d.mass-s.initialMass-s.budget;
  $('metrics-'+id).innerHTML=`<div><span>L₁ error</span><b>${fmt(d.l1)}</b></div><div><span>RMS error</span><b>${fmt(d.rms)}</b></div><div><span>min u / max u</span><b>${fmt(d.min)} / ${fmt(d.max)}</b></div><div><span>C + 2D · actual step</span><b>${fmt(s.load)}</b></div><div><span>Conservation residual</span><b>${fmt(residual)}</b></div><div><span>Total variation</span><b>${fmt(d.tv)}</b></div>`;
  if(!cache.has(s.time))cache.set(s.time,Array.from({length:401},(_,i)=>[i/400,B.exactPoint(p,i/400,s.time)]));
  const exact=cache.get(s.time),num=Array.from(s.u,(v,i)=>[(i+.5)*p.dx,v]);let lo=baselo,hi=basehi;
  if($('scale').value==='fit'){lo=Math.min(lo,...s.u.filter(Number.isFinite));hi=Math.max(hi,...s.u.filter(Number.isFinite));const pad=.06*(hi-lo);lo-=pad;hi+=pad;}
  const ch=chart('profile-'+id,[0,1],[lo,hi],$('scale').value==='shared'?'u · shared axes':'u · fitted axes');
  if(ch){if(showInitial)ch.line(init,'#959b8b',[3,4]);if(B.referenceAvailable(p,s.time))ch.line(exact,'#175842',[6,4]);ch.line(num,colors[id],[],true);}
  const emax=d.linf&&Number.isFinite(d.linf)?Math.max(1e-13,d.linf*1.12):.01,ec=chart('error-'+id,[0,1],[-emax,emax],'cell-average error');
  if(ec){ec.line([[0,0],[1,0]],'#a2aa9b',[3,4]);if(d.errors)ec.line(Array.from(d.errors,(v,i)=>[(i+.5)*p.dx,v]),colors[id]);else{ec.c.fillStyle='#6e796e';ec.c.textAlign='center';ec.c.fillText('Exact reference unavailable after tₛ',ec.X(.5),ec.Y(0)-8);}}
  $('note-'+id).textContent=!B.referenceAvailable(p,s.time)?'Past shock formation: the classical reference is unavailable; error norms are omitted.':d.min<lo||d.max>hi?'Profile exceeds shared axes. Use “Fit each method” to inspect its amplitude.':unresolved?'Viscous thickness < 3 cells: refine N before interpreting accuracy.':id==='maccormack'?'Bounded wiggles near a shock are not the same as growing instability.':`Errors compare cell averages. ${p.bc==='periodic'?'Mass should be constant.':'Mass changes through boundary fluxes.'}`;
 }
 if(n){const frozen=ids.filter(id=>states[id].frozen);$('run-message').className='run-message'+(frozen.length?' warn':'');$('run-message').textContent=n>=20000?'Step limit reached (20,000). Reset or use a shorter final time.':frozen.length?`${frozen.map(id=>B.methods[id]).join(', ')} frozen at their recorded times. Other methods continue.`:!B.referenceAvailable(p,t)?`Shock has formed. Numerical solutions continue; the classical exact reference and its error norms stop at tₛ = ${fmt(p.tShock)}.`:t>=end?'Final time reached. Compare errors, export the arrays, or refine the grid.':'Same grid, same time step. Inspect the actual C + 2D as the wave evolves.';}
 drawCharacteristics();drawHistory();
}
function drawHistory(){const vals=history.flatMap(z=>ids.map(id=>z[id])).filter(v=>v>0&&Number.isFinite(v));const lo=vals.length?Math.floor(Math.log10(Math.min(...vals)))-.5:-6,hi=vals.length?Math.max(lo+1,Math.ceil(Math.log10(Math.max(...vals)))):0;const ch=chart('history',[0,end],[lo,hi],'log₁₀ L₁ error vs physical time');if(ch)for(const id of ids)ch.line(history.map(z=>[z.time,z[id]>0?Math.log10(z[id]):null]),colors[id]);}
function drawCharacteristics(){
 const ch=chart('characteristics',[0,1],[0,end],'x horizontal · time t vertical');if(!ch)return;
 const paths=[];
 if(p.smooth){
  for(let j=0;j<=18;j++){const x0=j/18,v=B.initialPoint(p,x0);for(const shift of [-1,0,1])paths.push([[x0+shift,0],[x0+shift+v*end,end]]);}
  for(const a of paths)ch.line(a,'#9ab4a1');if(end>p.tShock)ch.line([[0,p.tShock],[1,p.tShock]],'#bf603b',[4,4]);
  $('wave-note').textContent=`Straight characteristics x = ξ + u₀(ξ)t first cross at tₛ = ${fmt(p.tShock)} (orange dashed). Their continuation beyond that time is multivalued, not the entropy solution. Neighboring periodic copies are shown.`;
 }else{
  const a=p.left,b=p.right,x0=p.x0;
  for(let j=0;j<=20;j++){const x=j/20,v=x<x0?a:b;ch.line([[x,0],[x+v*end,end]],'#a9bfaa');}
  if(a>b){ch.line([[x0,0],[x0+.5*(a+b)*end,end]],'#bf603b');$('wave-note').textContent=p.nu?`Inviscid characteristic guide (green); viscous shock center (orange). Speed s = ${(a+b)/2}; thickness scale 4ν/Δu = ${fmt(4*p.nu/(a-b))}.`:`Shock speed s = ${(a+b)/2}. Orange is xₛ(t) = ${x0} + st. Green paths are truncated conceptually at the shock; their straight continuations are only a guide.`;}
  else{for(let j=0;j<=10;j++)ch.line([[x0,0],[x0+(a+(b-a)*j/10)*end,end]],'#3b8e72');$('wave-note').textContent=`Fan edges: xL = ${x0} + (${a})t, xR = ${x0} + (${b})t. Inside the fan u = (x−x₀)/t. The entropy solution spreads.`;}
 }
 ch.line([[0,t],[1,t]],'#93657d',[3,4]);
}
function frame(now){if(!running)return;credit+=(now-lastFrame)/1000*+$('speed').value;lastFrame=now;const start=performance.now();while(credit>=1&&running&&performance.now()-start<22){if(!advance()){pause();break;}credit--;}credit=Math.min(credit,4);draw();if(running)requestAnimationFrame(frame);}
$('play').onclick=()=>{if(running){pause();return;}if(t>=end)reset();running=true;$('play').textContent='Pause';lastFrame=performance.now();credit=0;requestAnimationFrame(frame);};
$('step').onclick=()=>{pause();advance();draw();};$('reset').onclick=reset;
$('case').onchange=()=>{
 const meta=B.cases[$('case').value];if(meta.smooth)$('amplitude').value=meta.amplitude;
 $('end').value=B.problem($('case').value,+$('grid').value,+$('viscosity').value,+$('amplitude').value).end;reset();
};
$('amplitude').oninput=()=>{$('amplitude-value').textContent=$('amplitude').value;};
$('amplitude').onchange=()=>{$('end').value=B.problem($('case').value,+$('grid').value,+$('viscosity').value,+$('amplitude').value).end;reset();};
for(const id of ['grid','viscosity','cfl','timestep','end'])$(id).onchange=()=>{document.querySelectorAll('[data-preset]').forEach(b=>b.classList.remove('active'));reset();};
$('grid').oninput=()=>{$('grid-value').textContent=$('grid').value;};$('cfl').oninput=()=>{$('cfl-value').textContent=$('cfl').value;};
for(const id of ['scale','initial'])$(id).onchange=draw;
document.querySelectorAll('[data-preset]').forEach(button=>button.onclick=()=>{const key=button.dataset.preset;$('case').value=key==='unsafe'?'shock':key==='smooth'?'sine':key;$('amplitude').value=.5;$('end').value=key==='smooth'?.51:key==='entropy'?.3:.4;$('cfl').value=key==='unsafe'?1.3:.45;$('timestep').value=key==='unsafe'?'fixed':'adaptive';$('viscosity').value=.01;$('grid').value=key==='unsafe'?80:160;document.querySelectorAll('[data-preset]').forEach(b=>b.classList.toggle('active',b===button));reset();});
function csv(){const rows=[['method','time','step','x_left','x_right','cell_average','exact_cell_average','signed_error']];for(const id of ids){const s=states[id];for(let i=0;i<p.N;i++){const ref=B.exactAverage(p,i*p.dx,(i+1)*p.dx,s.time);rows.push([id,s.time,s.steps,i*p.dx,(i+1)*p.dx,s.u[i],ref??'',ref===null?'':s.u[i]-ref]);}}return rows.map(r=>r.join(',')).join('\n');}
$('download').onclick=()=>{const url=URL.createObjectURL(new Blob([csv()],{type:'text/csv'})),a=document.createElement('a');a.href=url;a.download=`burgers-${p.key}-t${t.toFixed(4)}.csv`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
$('primer-scheme').onchange=()=>{renderPrimerExample();$('primer-next').hidden=true;$('primer-step').textContent='Compute the next time level ↓';};
$('primer-step').onclick=()=>{renderPrimerExample();$('primer-next').hidden=!$('primer-next').hidden;$('primer-step').textContent=$('primer-next').hidden?'Compute the next time level ↓':'Hide the computed time level ↑';};
renderPrimerExample();new ResizeObserver(()=>{if(p)draw();}).observe($('plots'));reset();
