'use strict';
window.addEventListener('error',function(event){
 const box=document.getElementById('toast');
 if(box){box.hidden=false;box.textContent='화면 초기화 오류: '+(event.message||'화면 파일을 불러오지 못했습니다.')+' 새로고침 후 계속되면 실행 로그를 확인해 주세요.';}
});
