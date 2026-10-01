const node=(tag,text,className)=>{const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(className)n.className=className;return n;};
const link=(text,url)=>{const a=node('a',text);a.href=url;a.target='_blank';a.rel='noopener noreferrer';return a;};
const notice=text=>node('div',text,'notice');
async function api(path){const response=await fetch(path);const value=await response.json();if(!response.ok)throw Error(value.detail||'자료 조회 실패');return value;}
document.addEventListener('DOMContentLoaded',async()=>{
 const container=document.getElementById('saved-review-content');
 try{const id=location.pathname.split('/').pop();const dossier=await api('/api/dossiers/'+encodeURIComponent(id));container.replaceChildren();await renderSavedResults(dossier,container);}
 catch(e){container.replaceChildren(notice(e.message));}
});
