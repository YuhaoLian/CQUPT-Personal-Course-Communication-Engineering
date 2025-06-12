
    var hour=document.getElementById('h');//获取带有指定ID名称：h的元素
    var minute=document.getElementById('m');
    var second=document.getElementById('s');
    //获取当前时间
    function getCurrentTime(){
        var date=new Date();   //得到最新的时间
        var h=date.getHours();  //拿到小时数
        var m=date.getMinutes();    //拿到分钟数
        var s=date.getSeconds();    //拿到秒数
                
        if(h < 10) h='0'+h;//确保0-9时也显示成两位数
        if(m < 10) m='0'+m;//确保0-9分钟也显示成两位数
        if(s < 10) s='0'+s;//确保0-9秒也显示成两位数
                
        hour.innerHTML=h; //设置innerHTML指定的元素内容
        minute.innerHTML=m;
        second.innerHTML=s;
    }
    //每秒更新一次时间
    setInterval('getCurrentTime()',1000);        
