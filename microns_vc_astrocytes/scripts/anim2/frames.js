// node frames.js <outdir> <fps> <T> [i0 i1]  — renders anim.html frame by frame
const { chromium } = require('playwright');
const fs = require('fs');
(async () => {
  const [out, fps, T, i0, i1] = process.argv.slice(2);
  fs.mkdirSync(out, { recursive: true });
  const b = await chromium.launch({ executablePath: '/opt/pw-browsers/chromium',
    args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader'] });
  const p = await b.newPage({ viewport: { width: 1920, height: 1080 } });
  p.on('pageerror', e => console.log('PAGEERR', e.message));
  p.on('console', m => { if (m.type() === 'error') console.log('PAGE', m.text()); });
  await p.goto('http://localhost:8799/anim.html');
  await p.waitForFunction(() => window.READY === true, { timeout: 120000 });
  await p.waitForTimeout(1200);               // let webfonts settle
  const a = i0 === undefined ? 0 : Number(i0);
  const z = i1 === undefined ? Math.round(Number(T) * Number(fps)) : Number(i1);
  const t0 = Date.now();
  for (let k = a; k < z; k++) {
    await p.evaluate(t => window.renderAt(t), k / Number(fps));
    await p.screenshot({ path: `${out}/f${String(k).padStart(5, '0')}.png` });
    if ((k - a) % 60 === 0) console.log(k, z, ((Date.now() - t0) / 1000).toFixed(0) + 's');
  }
  await b.close();
})();
