var subm=1;
var loaded=0;

function openW (ERL,TARGET,OPT) 
{
  window.onerror=null;
  window.open(ERL,TARGET,OPT);
}

function combinedCheck()
{
    var one = Check();
    var two = Loaded();
    var combined = one && two;
    return  combined;
}

function Loaded()
{
    if (loaded)
        return true;
    else{
        if(subm)
            alert("Please wait until the entire page has been loaded.");
        return false;
    }
}       

function Check()
{
  if (subm)
     return true;
  else
     return false;
}

function SetChecks(f,select)
{
    
        for (var bo=2; bo<SetChecks.arguments.length;bo++)
        {
            var name = SetChecks.arguments[bo];
            for (var i=0;i<f.elements.length;i++)
            {
                var e=f.elements[i];
                if (e.name==name)
                {
                    if (select)
                    {
                        e.checked=true;
                    }
                    else
                    {
                        e.checked=false;
                    }
                }
            }
        }
    
    subm=0;
}


function ToggleChecks(f,hv,tn)
{
    var nv=true;
    for (var i=0;i<f.elements.length;i++)
    {
        var e=f.elements[i];
        if (e.name==hv)
        {
            if (e.value=="1")
            {
                e.value="0";
                nv=false;
            }
            else
            {
                e.value="1";
                nv=true;
            }
        }
    }
    for (var i=0;i<f.elements.length;i++)
    {
        var e=f.elements[i];
        if (e.name==tn)
        {
            e.checked=nv;
        }
    }
    subm=0;
}


function resizeFrame(fnum,nv)
{
    var gfs=parent.document.getElementById("gfs");
    var arr = new Array();
    arr = gfs.rows.split(',');
    arr[fnum*1] = nv;
    gfs.rows = arr[0]+','+arr[1]+','+arr[2];
}

function fullInput(f,obnm,msg)
{
    ret = 0;
    for (var i=0; i < f.elements.length; i++)
    {
        var e=f.elements[i];
        if (e.name == 'src2merge:list') {
            continue;
        }
        if ((obnm=='' && e.name.split(':lis').length==2) || e.name==obnm)
        {
            if (e.checked==true)
            {
                ret=1;   
            }
        }
    }
    if (ret==0) {
        alert(msg);
    }
    return ret;
}

function fullInput2(e,msg)
{    
    // see if e.value is not whitespaces only
    msg = "Keywords must be at least 3 characters long. Please provide appropriate keywords.";
    var nws = e.value.match(/\w+/);
    if (nws == null)
    {
        alert(msg);
    }
    else{
        var s = 1;
        var row_array = e.value.split('\n');
        for (var i=0;i<row_array.length;i=i+1){
           
            var row = row_array[i];
            var isnull = row.match(/\w+/);
            if (isnull != null) {
                var word_array = row.split(' ');
                s2 = 0;
                for(var j=0;j<word_array.length;j=j+1){
                    var word = word_array[j];
                    if(word.length >2){
                        s2 = 1;
                    }
                }
                if(s2 == 0){
                    s = 0;
                }
            }
        }
        if(s == 1){
            subm = 1;
        }
        else{
            alert(msg);
        }
    }
    
}

// Determine browser and version.

function Browser() {
  var ua, s, i;
  this.isIE    	=	false;
  this.isNS    	= 	false;
  this.isSafari =	false;
  this.isOpera  =   false;
  this.version 	= 	null;
  ua = navigator.userAgent;
  s = "MSIE";
  if ((i = ua.indexOf(s)) >= 0) {
    this.isIE = true;
    this.version = parseFloat(ua.substr(i + s.length));
    return;
  }
  s = "Opera";
  if ((i = ua.indexOf(s)) >= 0) {
    this.isOpera = true;
    return;
  }
  s = "Netscape6/";
  if ((i = ua.indexOf(s)) >= 0) {
    this.isNS = true;
    this.version = parseFloat(ua.substr(i + s.length));
    return;
  }
  s = "Safari";
  if ((i = ua.indexOf(s)) >= 0) {
    this.isSafari = true;
    this.version = parseFloat(ua.substr(i + s.length));
    return;
  }
  // Treat any other "Gecko" browser as NS 6.1.
  s = "Gecko";
  if ((i = ua.indexOf(s)) >= 0) {
    this.isNS = true;
    this.version = 6.1;
    return;
  }
}

var browser = new Browser();
//------------------------------------
// Move image with mousedrag functions
//------------------------------------

// Global object to hold scroll information.
var scrollObj = new Object();

function jsSetFocus(dvname,x,y) {
  var dv = document.getElementById(dvname);
  if (browser.isIE) {
    dv.scrollLeft = x - document.body.clientWidth/2;
    dv.scrollTop  = y - document.body.clientHeight/2;
  } else 
  {
    dv.scrollLeft = x - window.innerWidth/2;
    dv.scrollTop  = y - window.innerHeight/2;
  }
}

function jsSetFocus2(dvname,x,y) {
  var dv = parent.frames[1].document.getElementById(dvname);
  if (browser.isIE) {
    dv.scrollLeft = x - parent.frames[1].document.body.clientWidth/2;
    dv.scrollTop  = y - parent.frames[1].document.body.clientHeight/2;
  } else 
  {
    dv.scrollLeft = x - parent.frames[1].innerWidth/2;
    dv.scrollTop  = y - parent.frames[1].innerHeight/2;
  }
}

function scrollImage(event) {
  var el;
  var x, y;
  // Id (ImageLayer) of element that is going to be scrolled
  scrollObj.elNode = document.getElementById('imagediv');
  // Get cursor position with respect to the page.
  if (browser.isIE) {
    x = window.event.clientX + document.documentElement.scrollLeft
      + document.body.scrollLeft;
    y = window.event.clientY + document.documentElement.scrollTop
      + document.body.scrollTop;
  } else {
    x = event.clientX + window.scrollX;
    y = event.clientY + window.scrollY;
  }
  // Save starting positions of cursor and element.
  scrollObj.cursorStartX = x;
  scrollObj.cursorStartY = y;
  scrollObj.elStartLeft  = parseInt(scrollObj.elNode.scrollLeft, 10);
  scrollObj.elStartTop   = parseInt(scrollObj.elNode.scrollTop,  10);

  if (isNaN(scrollObj.elStartLeft)) scrollObj.elStartLeft = 0;
  if (isNaN(scrollObj.elStartTop))  scrollObj.elStartTop  = 0;

  // Capture mousemove and mouseup events on the page.
  if (browser.isIE) {
    document.attachEvent("onmousemove", dragStart);
    document.attachEvent("onmouseup",   dragStop);
    window.event.cancelBubble = true;
    window.event.returnValue = false;
  } else {
    document.addEventListener("mousemove", dragStart, 	true);
    document.addEventListener("mouseup",   dragStop, 	true);
    event.preventDefault();
  }
}


function dragStart(event) {
  var x, y, x_value, y_value;

  // The image
  var theImage		= document.getElementById('theimage');

  // Set cursor to grab (x-browser)
  theImage.className = 'xgrab';

  // Get cursor position with respect to the page.
  if (browser.isIE) {
    x = window.event.clientX + document.documentElement.scrollLeft
      + document.body.scrollLeft;
    y = window.event.clientY + document.documentElement.scrollTop
      + document.body.scrollTop;
  } else {
    x = event.clientX + window.scrollX;
    y = event.clientY + window.scrollY;
  }

  // Move scroll element by the same amount the cursor has moved.

  scrollObj.elNode.scrollLeft = (scrollObj.elStartLeft - x + scrollObj.cursorStartX);
  scrollObj.elNode.scrollTop  = (scrollObj.elStartTop  - y + scrollObj.cursorStartY);

  if (browser.isIE) {
    window.event.cancelBubble 	= true;
    window.event.returnValue 	= false;
  } else {
    event.preventDefault();
  }
}

function dragStop(event) {
  // The image
  var theImage		= document.getElementById('theimage');

  // Set cursor to openhand (x-browser)
  theImage.className = 'xopenhand';

  // Stop capturing mousemove and mouseup events.

  if (browser.isIE) {
    document.detachEvent("onmousemove", dragStart);
    document.detachEvent("onmouseup",   dragStop);
  } else {
    document.removeEventListener("mousemove", dragStart,   	true);
    document.removeEventListener("mouseup",   dragStop, 	true);
  }
}

function TrimString(sInString) {
  //sInString = sInString.replace( /^\s+/g, "" );
  sInString = sInString.replace(/\s+/g, "" )
  return sInString;}


function constraintInputCheck() {
 for (var r=1; r<3; r++){
  if(document.getElementById("nn_"+r).checked){
   mr = document.getElementById("mems_"+r).value;
   m = parseInt(mr);
   if (isNaN(m) || m < 2) {
    document.getElementById("mems_"+r).value = "2";
    document.getElementById("mems_"+r).style.background="#ffaaaa";
    document.getElementById("mems_"+r).style.color="#000000";
   }
   else if (m+"" != mr) {
    document.getElementById("mems_"+r).value = m+"";
    document.getElementById("mems_"+r).style.background="#ffaaaa";
    document.getElementById("mems_"+r).style.color="#000000";             
   }

   cr = document.getElementById("cc_"+r).value;
   c = parseFloat(cr);
   if (isNaN(c) || c < 0 || c > 1) {
    document.getElementById("cc_"+r).value = "0.0";
    document.getElementById("cc_"+r).style.background="#ffaaaa";
    document.getElementById("cc_"+r).style.color="#000000";
   }
   else if (cr != "0.0" && c+"" != cr) {
    document.getElementById("cc_"+r).value = c+"";
    document.getElementById("cc_"+r).style.background="#ffaaaa";
    document.getElementById("cc_"+r).style.color="#000000";
   }

   mc = document.getElementById("mc_"+r).value;
   m = parseInt(mc);
   if (isNaN(m) || m < 1) {
    document.getElementById("mc_"+r).value = "1";
    document.getElementById("mc_"+r).style.background="#ffaaaa";
    document.getElementById("mc_"+r).style.color="#000000";
   }
   else if (m+"" != mc) {
    document.getElementById("mc_"+r).value = m+"";
    document.getElementById("mc_"+r).style.background="#ffaaaa";
    document.getElementById("mc_"+r).style.color="#000000";
   }

   pc = document.getElementById("pcut_"+r).value;
   p = parseFloat(pc);
   if (isNaN(p) || p < 0 || p > 1) {
    document.getElementById("pcut_"+r).value = "0.001";
    document.getElementById("pcut_"+r).style.background="#ffaaaa";
    document.getElementById("pcut_"+r).style.color="#000000";
   }
   else if (pc != "0.001" && p+"" != pc) {
    document.getElementById("pcut_"+r).value = p+"";
    document.getElementById("pcut_"+r).style.background="#ffaaaa";
    document.getElementById("pcut_"+r).style.color="#000000";
   }
  }
 }

 pc = document.getElementById("pcut_pwy").value;
 p = parseFloat(pc);
 if (isNaN(p) || p < 0 || p > 1) {
  document.getElementById("pcut_pwy").value = "0.01";
  document.getElementById("pcut_pwy").style.background="#ffaaaa";
  document.getElementById("pcut_pwy").style.color="#000000";
 }
 else if (pc != "0.01" && p+"" != pc) {
  document.getElementById("pcut_pwy").value = p+"";
  document.getElementById("pcut_pwy").style.background="#ffaaaa";
  document.getElementById("pcut_pwy").style.color="#000000";
 }

 mc = document.getElementById("mc_p").value;
 m = parseInt(mc);
 if (isNaN(m) || m < 1) {
  document.getElementById("mc_p").value = "1";
  document.getElementById("mc_p").style.background="#ffaaaa";
  document.getElementById("mc_p").style.color="#000000";
 }
 else if (m+"" != mc) {
  document.getElementById("mc_p").value = m+"";
  document.getElementById("mc_p").style.background="#ffaaaa";
  document.getElementById("mc_p").style.color="#000000";
 }

 go2 = document.getElementById("pcut_go2").value;
 p = parseFloat(go2);
 if (isNaN(p) || p < 0 || p > 1) {
  document.getElementById("pcut_go2").value = "0.01";
  document.getElementById("pcut_go2").style.background="#ffaaaa";
  document.getElementById("pcut_go2").style.color="#000000";
 }
 else if (go2 != "0.01" && p+"" != go2) {
  document.getElementById("pcut_go2").value = p+"";
  document.getElementById("pcut_go2").style.background="#ffaaaa";
  document.getElementById("pcut_go2").style.color="#000000";
 }

 go3 = document.getElementById("pcut_go3").value;
 p = parseFloat(go3);
 if (isNaN(p) || p < 0 || p > 1) {
  document.getElementById("pcut_go3").value = "0.01";
  document.getElementById("pcut_go3").style.background="#ffaaaa";
  document.getElementById("pcut_go3").style.color="#000000";
 }
 else if (go3 != "0.01" && p+"" != go3) {
  document.getElementById("pcut_go3").value = p+"";
  document.getElementById("pcut_go3").style.background="#ffaaaa";
  document.getElementById("pcut_go3").style.color="#000000";
 }
}

getVizCookie = function (){
    var a = document.cookie;
    var v = -1;
    if(a){
        var lines = a.split(";");
        for (var i in lines) {
            kv = lines[i].replace(/^\s+|\s+$/g,"").split("=");
            if (kv[0] == "vizCookie"){
                v = kv[1];
                break;
            }
        }
    }
    if (v == -1) {
        v = "java";
    }
    // browser check; change if IE
    if (browser.isIE) {
        v = "python";
    }
    document.getElementById("viztype_"+v).checked = true;
}


function setVizCookie(vistype){
    var a = new Date();
    a = new Date(a.getTime() + 1000*60*60*24*365);
    document.cookie = 'vizCookie='+vistype+'; expires='+a.toGMTString()+';';
    //alert(document.cookie);
}

function highlightCell(rowids){
    var clicked = document.getElementById('clicked').value;
    if(clicked != ""){
        var highlightedRows = clicked.split(';;');
        //alert(highlightedRows);
        //alert(highlightedRows.length);
        for(var i=0;i<highlightedRows.length;i=i+1){
           document.getElementById(highlightedRows[i]).style.color='black';
        }
    }
    var rows = rowids.split(';;');
    for(var i=0;i<rows.length;i=i+1){
        document.getElementById(rows[i]).style.color='red';
    }
    document.getElementById('clicked').value=rowids;
}

function selectEsets(tp,b) {
    var cb = document.getElementsByName('objs:list');
    for (var c in cb) {
        c = cb[c]
        if (c.value[0] == tp) {
            c.checked = b;
        }
    }
}
