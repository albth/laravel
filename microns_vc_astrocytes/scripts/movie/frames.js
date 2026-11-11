// usage: node frames.js <L> <R> <cfg.json> <outdir> <fps> <T> [only_times_csv]
const { chromium } = require('playwright'); const fs=require('fs');
(async()=>{const [L,R,cfgf,out,fps,T,only]=process.argv.slice(2);
 const cfg=JSON.parse(fs.readFileSync(cfgf,'utf8')); fs.mkdirSync(out,{recursive:true});
 const b=await chromium.launch({executablePath:'/opt/pw-browsers/chromium',args:['--use-angle=swiftshader','--enable-unsafe-swiftshader','--ignore-gpu-blocklist']});
 const p=await b.newPage({viewport:{width:1920,height:1080}});
 p.on('console',m=>{if(m.type()==='error')console.log('PAGE',m.text())}); p.on('pageerror',e=>console.log('PAGEERR',e.message));
 await p.goto(`http://localhost:8765/movie.html?L=${L}&R=${R}`); await p.waitForFunction(()=>window.READY===true,{timeout:300000});
 await p.evaluate(c=>window.setup(c),cfg);
 const times= only? only.split(',').map(Number) : [...Array(Math.round(T*fps)).keys()].map(i=>i/fps);
 let k=0; const t0=Date.now();
 for(const t of times){ await p.evaluate(t=>window.renderAt(t),t); await p.screenshot({path:`${out}/f${String(k).padStart(5,'0')}.png`}); k++;
   if(k%60==0) console.log(k,times.length,((Date.now()-t0)/1000).toFixed(0)+'s'); }
 await b.close();})();
