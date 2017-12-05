(function() {
    var requestAnimationFrame = window.requestAnimationFrame || window.mozRequestAnimationFrame || window.webkitRequestAnimationFrame || window.msRequestAnimationFrame;
    window.requestAnimationFrame = requestAnimationFrame;
})();
//<canvas id="glcanvas" width="10" height="10" style="border:10px solid #000000;"></canvas>
var n = 0, a = 0, b=0, c=0, d=0, e=0, f=0, g=0, i=0;

hun = document.getElementById("hundred");
var canvas = document.getElementById("text");
var ctx=canvas.getContext("2d");
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;
  camera = new THREE.PerspectiveCamera( 60, window.innerWidth / window.innerHeight, 1, 30000 );
  controls = new THREE.FirstPersonControls( camera );
w=canvas.width;
h=canvas.height;
ctx.clearRect(0,0,w,h);
//ctx.setLineDash([5, 0]);
keys = [];
var y1=y2=0.8*h,y3=0.957*h,y4=0.8857*h;
var nei=controls.a(), stei=controls.b();

function update(){
  var c1="white",c2="white",xs=0.10*w,ys=0.7857*h;
  ctx.clearRect(0,0,w,h);
  ctx.lineWidth = 0.001*w;
  ctx.strokeStyle="black";

  controls.schub = (0.8*h-y1)/0.15/h;
  controls.update( clock.getDelta() );

//Cockpit-Rahmen
  ctx.beginPath();
  ctx.moveTo(0,0);
  ctx.lineTo(w,0);
  ctx.lineTo(w,h);
  ctx.lineTo(0.974*w,h);
  ctx.lineTo(0.974*w,0.67*h);
  ctx.lineTo(0.93*w,0.561*h);
  ctx.lineTo(0.9125*w,0.561*h);
  ctx.lineTo(0.94*w,0.026*h);
  ctx.lineTo(0.508*w,0.066*h);
  ctx.lineTo(0.512*w,0.49*h);
  ctx.lineTo(0.546*w,0.561*h);
  ctx.lineTo(0.454*w,0.561*h);
  ctx.lineTo(0.488*w,0.49*h);
  ctx.lineTo(0.492*w,0.066*h);
  ctx.lineTo(0.06*w,0.026*h);
  ctx.lineTo(0.0875*w,0.561*h);
  ctx.lineTo(0.07*w,0.561*h);
  ctx.lineTo(0.026*w,0.67*h);
  ctx.lineTo(0.026*w,h);
  ctx.lineTo(0,h);
  ctx.lineTo(0,0);
  ctx.fillStyle=" #ccccb3";
  ctx.fill();
  ctx.closePath();

//Untere Box für Instrumente
  ctx.beginPath();
  ctx.moveTo(0.026*w,0.67*h);
  ctx.lineTo(0.974*w,0.67*h);
  ctx.lineTo(0.974*w,h);
  ctx.lineTo(0.026*w,h);
  ctx.lineTo(0.026*w,0.67*h);
  ctx.fillStyle=" #A4A4A4";
  ctx.fill();
  ctx.closePath();
//obere Box für Knöpfe

  ctx.beginPath();
  ctx.moveTo(0.07*w,0.561*h);
  ctx.lineTo(0.93*w,0.561*h);
  ctx.lineTo(0.974*w,0.67*h);
  ctx.lineTo(0.026*w,0.67*h);
  ctx.lineTo(0.07*w,0.561*h);


  ctx.closePath();
  ctx.fill();
  //Höhenmesser
  fla();

/*
  showN();
  showN3();
  showN5();
  showN7();
  showN9();
  showN11();
  showN13();
  showN14();
  showN15();
  */
            if (keys[39]) {
       // right arrow
       c2="red";
       xs=0.13*w;
   }
   if (keys[37]) {
        // left arrow
  c1="red";
  xs=0.1*w;
   }
   if (keys[40]) {
// up arrow
ys=0.807*h;
}
if (keys[38]) {
// down arrow
ys=0.764*h;
}


ctx.beginPath();
ctx.lineWidth="2";
ctx.strokeStyle = "white";
//ctx.setLineDash([5, 0]);
/*Schubkontrolle*/
ctx.strokeRect(0.46*w,0.714*h,0.03*w,0.171*h);
ctx.strokeRect(0.51*w,0.714*h,0.03*w,0.171*h);
/*Landeklappenhebel+Bremsklappenhebel*/
ctx.strokeRect(0.46*w,0.929*h,0.02*w,0.057*h);
ctx.strokeRect(0.51*w,0.929*h,0.02*w,0.057*h);
/*Gearlevel+Warnanzeige*/
ctx.strokeRect(0.82*w,0.857*h,0.026*w,0.1*h);
ctx.strokeRect(0.88*w,0.857*h,0.07*w,0.1*h);

ctx.fillStyle= "black";
/*Sidestick*/
ctx.arc(xs,ys,0.0571*h,0,2*Math.PI);
/*Pedale*/
ctx.rect(0.11*w,0.9*h,0.04*w,0.0714*h);
ctx.arc(0.15*w,0.9357*h,0.0357*h,0,2*Math.PI);
ctx.rect(0.06*w,0.9*h,0.04*w,0.0714*h);
ctx.arc(0.06*w,0.9357*h,0.0357*h,0,2*Math.PI);
/*Schubkontrolle*/
ctx.rect(0.455*w,y1,0.04*w,0.02857*h);
ctx.rect(0.505*w,y2,0.04*w,0.02857*h);
/*Landeklappenhebel+Bremsklappenhebel*/
ctx.rect(0.455*w,y3,0.030*w,0.0214*h);
ctx.arc(0.52*w,0.957*h,0.0214*h,0,2*Math.PI);
/*Gearlevel+Warnanzeige*/
ctx.rect(0.815*w,y4,0.036*w,0.0214*h);
ctx.fill();



drawArrow(0.12*w,0.9357*h,0.15*w,0.9357*h,c2);
drawArrow(0.09*w,0.9357*h,0.06*w,0.9357*h,c1);
ctx.strokeStyle = "white";

ctx.closePath();
/*
//Koordinatenanzeige
ctx.fillStyle = "grey";
ctx.fillRect(0.7*w,0.7*h,0.1*w,0.1*w);
ctx.fillStyle="black";
ctx.fillText("Position : ", 0.75*w, 0.72*h);
ctx.fillText("Z : ",0.75*w,0.75*h);
ctx.fillText(Math.round(camera.position.z), 0.79*w, 0.75*h);
ctx.fillText("X : ",0.75*w,0.78*h);
ctx.fillText(Math.round(camera.position.x), 0.79*w, 0.78*h);
ctx.fillText("Y : ",0.75*w,0.81*h);
ctx.textAlign = "right";
ctx.fillText(Math.round(camera.position.y), 0.79*w, 0.81*h);
*/
   requestAnimationFrame(update);
}

var click1,click2,click3,click4,mx,my;

document.body.addEventListener("keydown", function(e) {
    keys[e.keyCode] = true;
});

document.body.addEventListener("keyup", function(e) {
    keys[e.keyCode] = false;
});
document.body.addEventListener("mousemove", function(e) {
if(click1){
  if(my>0.7285*h&&my<0.8714*h){
    y1=y1+e.offsetY-my;
  }
}
  if(click2){
    if(my>0.7285*h&&my<0.8714*h){
      y2=y2+e.offsetY-my;
    }
  }
  if(click3){
    if(my>0.9357*h&&my<0.9785*h){
      y3=y3+e.offsetY-my;
    }
  }
  if(click4){
    if(my>0.8642*h&&my<0.95*h){
      y4=y4+e.offsetY-my;
    }
  }
  my=e.offsetY;

});
document.body.addEventListener("mousedown", function(e) {
  mx=e.clientX;
  my=e.offsetY;
  if(mx > 0.455*w && mx < 0.495*w && my>y1 && my<(y1+0.0286*w)){
    click1=true;
  }
  if(mx > 0.505*w && mx < 0.545*w && my>y2 && my<(y2+0.0286*w)){
    click2=true;
  }
  if(mx > 0.435*w && mx < 0.465*w && my>y3 && my<(y3+0.0214*w)){
    click3=true;
  }
  if(mx > 0.815*w && mx < 0.851*w && my>y4 && my<(y4+0.0214*w)){
    click4=true;
  }

});
document.body.addEventListener("mouseup", function(e) {
  click1=click2=click3=click4=false;
});

window.addEventListener("load",function(){
    update();
});

function drawArrow(fromx, fromy, tox, toy,color){
                //variables to be used when creating the arrow
                //ctx.setLineDash([5, 0]);
                var headlen = 1;

                var angle = Math.atan2(toy-fromy,tox-fromx);

                //starting path of the arrow from the start square to the end square and drawing the stroke
                ctx.beginPath();
                ctx.moveTo(fromx, fromy);
                ctx.lineTo(tox, toy);
                ctx.strokeStyle = color;
                ctx.lineWidth = 0.01*w;
                ctx.stroke();

                //starting a new path from the head of the arrow to one of the sides of the point
                ctx.beginPath();
                ctx.moveTo(tox, toy);
                ctx.lineTo(tox-headlen*Math.cos(angle-Math.PI/7),toy-headlen*Math.sin(angle-Math.PI/7));

                //path from the side point of the arrow, to the other side point
                ctx.lineTo(tox-headlen*Math.cos(angle+Math.PI/7),toy-headlen*Math.sin(angle+Math.PI/7));

                //path from the side point back to the tip of the arrow, and then again to the opposite side point
                ctx.lineTo(tox, toy);
                ctx.lineTo(tox-headlen*Math.cos(angle-Math.PI/7),toy-headlen*Math.sin(angle-Math.PI/7));

                //draws the paths created above
                ctx.strokeStyle = color;
                ctx.lineWidth = 0.01*w;
                ctx.stroke();
                ctx.fillStyle = color;
                ctx.fill();
            }

            function resizeCanvas() {
              canvas.width = window.innerWidth;
              canvas.height = window.innerHeight;
              w=canvas.width;
              h=canvas.height;
                // in this case just render when the window is resized.
                update();

            }

            window.addEventListener('resize', resizeCanvas, false);

            /*Knöpfe*/
            function showN() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "green";
                ctx.beginPath();
                ctx.arc(0.18*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }
            function showN2() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "red";
                ctx.beginPath();
                ctx.arc(0.18*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }
            function showN3() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "green";
                ctx.beginPath();
                ctx.arc(0.28*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }
            function showN4() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "red";
                ctx.beginPath();
                ctx.arc(0.28*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }
            function showN5() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "green";
                ctx.beginPath();
                ctx.arc(0.38*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }
            function showN6() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "red";
                ctx.beginPath();
                ctx.arc(0.38*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }
            function showN7() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "green";
                ctx.beginPath();
                ctx.arc(0.48*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }
            function showN8() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "red";
                ctx.beginPath();
                ctx.arc(0.48*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }
            function showN9() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "green";
                ctx.beginPath();
                ctx.arc(0.62*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }
            function showN10() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "red";
                ctx.beginPath();
                ctx.arc(0.62*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }
            function showN11() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "green";
                ctx.beginPath();
                ctx.arc(0.76*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }
            function showN12() {
                //ctx.clearRect(0, 0, can.width, can.height);
                ctx.fillStyle = "red";
                ctx.beginPath();
                ctx.arc(0.76*w,41/70*h,1/70*h,0,Math.PI * 2, true);
                ctx.closePath();
                ctx.fill();
                ctx.stroke();
            }


            function incr() {
                n++;
                if(n%2==0){
                  showN();
            }
               else { showN2();
               }
             }
             function incr2() {
                 a++;
                 if(a%2==0){
                   showN3();
             }
                else { showN4();
                }
              }
              function incr3() {
                  b++;
                  if(b%2==0){
                    showN5();
              }
                 else { showN6();
                 }
               }
               function incr4() {
                   c++;
                   if(c%2==0){
                     showN7();
               }
                  else { showN8();
                  }
                }
                function incr5() {
                    d++;
                    if(d%2==0){
                      showN9();
                }
                   else { showN10();
                   }
                 }
                 function incr6() {
                     e++;
                     if(e%2==0){
                       showN11();
                 }
                    else { showN12();
                    }
                  }


              function showN13() {

                  // large, centered, bright green text
                  ctx.font = "18pt Helvetica";
                  ctx.textAlign = "right";
                  ctx.textBaseline = "right";
                  ctx.fillStyle = "rgb(255,222,173)";
                  ctx.clearRect(0.5*w,0.571*h,0.08*w,1/7*h);
                  // draw text at center, max length to fit on canvas
                  ctx.fillText(f, 0.57*w, 0.6*h, w - 2);

              }
              function showN14() {

                  // large, centered, bright green text
                  ctx.font = "18pt Helvetica";
                  ctx.textAlign = "right";
                  ctx.textBaseline = "right";
                  ctx.fillStyle = "rgb(255,222,173)";
                  ctx.clearRect(0.64*w,4/7*h,0.08*w,1/7*h);
                  // draw text at center, max length to fit on canvas
                  ctx.fillText(g, 0.71*w, 0.6*h, w - 2);

              }
              function showN15() {

                  // large, centered, bright green text
                  ctx.font = "18pt Helvetica";
                  ctx.textAlign = "right";
                  ctx.textBaseline = "right";
                  ctx.fillStyle = "rgb(255,222,173)";
                  ctx.clearRect(0.78*w,4/7*h,0.08*w,1/7*h);
                  // draw text at center, max length to fit on canvas
                  ctx.fillText(i, 0.85*w, 42/70*h, w - 2);

              }

              function incr7() {
                f=f+100;
                showN13();
              }

             function decr1() {
               f=f-100;
               showN13();
              }
              function incr8() {
                g=g+10;
                showN14();
              }

             function decr2() {
               g=g-10;
               showN14();
              }
              function incr9() {
                i=i+10;
                showN15();
              }

             function decr3() {
               i=i-10;
               showN15();
              }

              function fla(){
              //Geschwindigkeit in Knoten (-> Verschiebung in y-Richtung: -250+4*v)
              var v = 180;
              //Höhe in Fuß
              nei=controls.a(); stei=0.005*h*controls.b();
              var a=h/10;
              ctx.translate(0.31*w,0.83*h);
              ctx.fillStyle = "lightblue"
              ctx.fillRect(-0.09*w,-0.15*h,0.18*w,0.3*h);
              ctx.translate(0,stei);
              ctx.rotate(-nei*Math.PI/180);

              ctx.fillStyle="lightgreen"
              ctx.fillRect(-0.09*w,0.0*h,0.18*w,0.3*h);
// Beschriftung der Skala vom fla

              ctx.font="9px Arial";
              ctx.fillStyle = "black";
              ctx.textAlign = "left";
              ctx.fillText("10 ", -0.039*w, -0.045*h);
              ctx.fillText("10", 0.031*w, -0.045*h);
              ctx.fillText("20 ", -0.039*w, -0.095*h);
              ctx.fillText("20", 0.031*w, -0.095*h);
              ctx.fillText("10 ", -0.039*w, 0.055*h);
              ctx.fillText("10", 0.031*w, 0.055*h);
              ctx.fillText("20 ", -0.039*w, 0.105*h);
              ctx.fillText("20", 0.031*w, 0.105*h);

              ctx.lineWidth = 0.0005*w;
              ctx.strokeStyle = "black";
// skala vom fla
              ctx.moveTo(-0.01*w,-0.0125*h);
              ctx.lineTo(0.01*w,-0.0125*h);
              ctx.moveTo(-0.02*w,-0.025*h);
              ctx.lineTo(0.02*w,-0.025*h);
              ctx.moveTo(-0.01*w,-0.0375*h);
              ctx.lineTo(0.01*w,-0.0375*h);
              ctx.moveTo(-0.03*w,-0.05*h);
              ctx.lineTo(0.03*w,-0.05*h);
              ctx.moveTo(-0.01*w,-0.0625*h);
              ctx.lineTo(0.01*w,-0.0625*h);
              ctx.moveTo(-0.02*w,-0.075*h);
              ctx.lineTo(0.02*w,-0.075*h);
              ctx.moveTo(-0.01*w,-0.0875*h);
              ctx.lineTo(0.01*w,-0.0875*h);
              ctx.moveTo(-0.03*w,-0.1*h);
              ctx.lineTo(0.03*w,-0.1*h);
              ctx.moveTo(-0.09*w,0);
              ctx.lineTo(0.09*w,0);

              ctx.moveTo(-0.01*w,0.0125*h);
              ctx.lineTo(0.01*w,0.0125*h);
              ctx.moveTo(-0.02*w,0.025*h);
              ctx.lineTo(0.02*w,0.025*h);
              ctx.moveTo(-0.01*w,0.0375*h);
              ctx.lineTo(0.01*w,0.0375*h);
              ctx.moveTo(-0.03*w,0.05*h);
              ctx.lineTo(0.03*w,0.05*h);
              ctx.moveTo(-0.01*w,0.0625*h);
              ctx.lineTo(0.01*w,0.0625*h);
              ctx.moveTo(-0.02*w,0.075*h);
              ctx.lineTo(0.02*w,0.075*h);
              ctx.moveTo(-0.01*w,0.0875*h);
              ctx.lineTo(0.01*w,0.0875*h);
              ctx.moveTo(-0.03*w,0.1*h);
              ctx.lineTo(0.03*w,0.1*h);
              ctx.stroke();




/*
              ctx.moveTo(-0.1188*w,0);
              ctx.lineTo(-0.11664*w,-0.141*h);
              ctx.moveTo(210,0);
              ctx.lineTo(213,15);
              ctx.moveTo(340,12);
              ctx.lineTo(333,28);
              ctx.moveTo(160,12);
              ctx.lineTo(167,28);
              ctx.moveTo(400,18);
              ctx.lineTo(375,55);
              ctx.moveTo(100,18);
              ctx.lineTo(125,55);
              ctx.moveTo(435,75);
              ctx.lineTo(424,85);
              ctx.moveTo(65,75);
              ctx.lineTo(76,85);
              ctx.stroke();
              */
              ctx.rotate(nei*Math.PI/180);
              ctx.translate(-0.31*w,-0.83*h-stei);
// konstante Anzeigen im fla
              ctx.translate(0.31*w,0.83*h);
              ctx.fillStyle="black";
              ctx.beginPath();
              ctx.moveTo(-0.005*w,-0.15*h);
              ctx.lineTo(0.005*w,-0.15*h);
              ctx.lineTo(0,-0.14*h);
              ctx.lineTo(-0.005*w,-0.15*h);
              ctx.fill();
              //dreieck für Querneigungswinkel
              ctx.moveTo(0,-0.14*h);
              ctx.lineTo(0.005*w,-0.13*h);
              ctx.lineTo(-0.005*w,-0.13*h);
              ctx.lineTo(0,-0.14*h);
              ctx.stroke();



              ctx.fillStyle="black";
              ctx.beginPath();
              ctx.moveTo(0.03*w,-0.005*h);
              ctx.lineTo(0.07*w,-0.005*h);
              ctx.lineTo(0.07*w,0.005*h);
              ctx.lineTo(0.037*w,0.005*h);
              ctx.lineTo(0.037*w,0.025*h);
              ctx.lineTo(0.03*w,0.025*h);
              ctx.lineTo(0.03*w,-0.005*h);
              ctx.moveTo(-0.03*w,-0.005*h);
              ctx.lineTo(-0.07*w,-0.005*h);
              ctx.lineTo(-0.07*w,0.005*h);
              ctx.lineTo(-0.037*w,0.005*h);
              ctx.lineTo(-0.037*w,0.025*h);
              ctx.lineTo(-0.03*w,0.025*h);
              ctx.lineTo(-0.03*w,-0.005*h);
              ctx.fill();

              ctx.fillStyle = "black"
              ctx.beginPath();
              ctx.arc(0,0,0.005*w,0,Math.PI * 2, true);
              ctx.closePath();
              ctx.fill();
              //ctx.translate(-0.31*w,-0.83*h);

//Bänder für Höhe und Geschwindigkeit
              ctx.fillStyle="#585858";
              ctx.fillRect(-0.13*w,-0.2*h,0.04*w,0.8*h);
              ctx.fillRect(0.09*w,-0.2*h,0.04*w,0.8*h);
              var Geschwindigkeit = controls.movementSpeed;

              ctx.translate(0,-490+2*(Math.round(Geschwindigkeit)));
              ctx.translate(21,0);
              ctx.font="12px Arial white";
              ctx.fillStyle = "white";
              ctx.textAlign = "right";
              var grenze1=Math.round(Math.round(Geschwindigkeit)/20)*2+6;
              for (var i=grenze1; i>-1;){
                //ctx.moveTo(-145,490-i*20);
                //ctx.lineTo(-130,490-i*20);
                //ctx.moveTo(-145,490-i*10);
                //ctx.lineTo(-130,490-i*10);
                ctx.fillText(i*10,-150,495-i*20);
                ctx.stroke();
                i=i-2;
              }
              for (var i=grenze1*2; i>-1;){
                ctx.moveTo(-145,490-i*10);
                ctx.lineTo(-130,490-i*10);
                i=i-2;
              }
              ctx.stroke();
              ctx.translate(-21,0);

              ctx.translate(0,+490-2*(Math.round(Geschwindigkeit)));


              ctx.translate(0,-490+2*(Math.round(camera.position.y)/10)); //h=0 bei translate(0,-490)
              ctx.translate(-21,0);
              ctx.font="12px Arial white";
              ctx.fillStyle = "white";
              ctx.textAlign = "left";
              var grenze2=Math.round(Math.round(camera.position.y)/200)*2+8;
              for (var i=grenze2; i>-1;){
              //  ctx.moveTo(145,490-i*20);
              //  ctx.lineTo(130,490-i*20);
              //  ctx.moveTo(145,490-i*10);
                //ctx.lineTo(130,490-i*10);
                ctx.fillText(i*100,150,495-i*20);
                i=i-2;
              }
              for (var i=grenze2*2; i>-1;){
                ctx.moveTo(145,490-i*10);
                ctx.lineTo(130,490-i*10);
                i=i-2;
              }
              ctx.stroke();
              ctx.translate(21,0);
ctx.translate(0,490-2*(Math.round(camera.position.y)/10));
              ctx.fillStyle="black";



ctx.translate(31,0);
//Anzeigeboxen für Geschwindigkeit und Höhe
              ctx.fill();
              ctx.fillStyle="white";
              ctx.beginPath();
              ctx.moveTo(-187,-12);
              ctx.lineTo(-149,-12);
              ctx.lineTo(-149,-7);
              ctx.lineTo(-144,1);
              ctx.lineTo(-149,9);
              ctx.lineTo(-149,14);
              ctx.lineTo(-187,14);
              ctx.closePath();
              ctx.fill();
ctx.translate(-31,0);
ctx.translate(-31,0);
              ctx.beginPath();
              ctx.moveTo(187,-12);
              ctx.lineTo(149,-12);
              ctx.lineTo(149,-7);
              ctx.lineTo(144,1);
              ctx.lineTo(149,9);
              ctx.lineTo(149,14);
              ctx.lineTo(187,14);
              ctx.closePath();
              ctx.fill();
ctx.translate(31,0);
              ctx.fillStyle="black";
ctx.translate(31,0);
              ctx.beginPath();
              ctx.moveTo(-187,-10);
              ctx.lineTo(-151,-10);
              ctx.lineTo(-151,-5);
              ctx.lineTo(-146,1);
              ctx.lineTo(-151,7);
              ctx.lineTo(-151,12);
              ctx.lineTo(-187,12);
              ctx.lineTo(-187,-10);
              ctx.closePath();
              ctx.fill();
ctx.translate(-31,0);

ctx.translate(-31,0);
              ctx.beginPath();
              ctx.moveTo(187,-10);
              ctx.lineTo(151,-10);
              ctx.lineTo(151,-5);
              ctx.lineTo(146,1);
              ctx.lineTo(151,7);
              ctx.lineTo(151,12);
              ctx.lineTo(187,12);
              ctx.closePath();
              ctx.fill();
              ctx.stroke();
              ctx.fillStyle="white"
ctx.translate(31,0);


ctx.translate(31,0);
              ctx.font="12px Arial white";
              ctx.fillStyle = "white";
              ctx.textAlign = "left";
              ctx.fillText(Math.round(Geschwindigkeit),-184,4);
              ctx.translate(-31,0);
              ctx.translate(-31,0);
              ctx.font="12px Arial white";
              ctx.fillStyle = "white";
              ctx.textAlign = "right";

              ctx.fillText(Math.round(camera.position.y),184,0.005*h);
              ctx.translate(31,0);
                ctx.translate(-0.31*w,-0.83*h);
              /*Umrandung der Höhenanzeige*/
              ctx.translate(0.31*w,0.83*h);
              ctx.fillStyle = "#A4A4A4"
              ctx.fillRect(-0.2*w,0.15*h,0.5*w,0.3*h);
              ctx.fillRect(-0.2*w,-0.25*h,0.5*w,0.1*h);
              //ctx.fillRect(-0.2*w,-0.15*h,0.11*w,0.3*h);
              //ctx.fillRect(0.09*w,-0.15*h,0.11*w,0.3*h);

              ctx.moveTo(-0.2*w,-0.16*h);
              ctx.lineTo(0.3*w,-0.16*h);
              ctx.stroke();
              ctx.translate(-0.31*w,-0.83*h);

//navigation
ctx.fillStyle="black";
ctx.fillRect(693,507,230,230);
//Ausrichtung
var a=0;
//heading
var b=90;
ctx.lineWidth=0.1;
ctx.translate(810,620);

ctx.rotate(-a*Math.PI/180);
ctx.beginPath();
ctx.arc(0,0,90,0,Math.PI * 2, true);
ctx.closePath();
ctx.strokeStyle="white";
ctx.stroke();
ctx.font="12px Arial";
ctx.fillStyle = "white";
ctx.textAlign = "left";
for(var i=0;i<36;){
ctx.moveTo(90,0);
ctx.lineTo(100,0);
ctx.strokeStyle="white";
ctx.stroke();
if(i%3==0){
  ctx.fillText(i,103,5);
}
ctx.rotate(-10*Math.PI/180);
i=i+1;
}
for(var n=1;n<74;){
ctx.rotate(-5*Math.PI/180);
ctx.moveTo(90,0);
ctx.lineTo(95,0);
ctx.strokeStyle="white";
ctx.stroke();
n=n+1;
}

ctx.lineWidth=1;
//ctx.setLineDash([5, 15]);
ctx.strokeStyle="white";
ctx.beginPath();
ctx.arc(0,0,60,0,Math.PI * 2, true);
ctx.closePath();
ctx.strokeStyle="white";
ctx.stroke();

ctx.rotate(a*Math.PI/180);
ctx.rotate(36*10*Math.PI/180);
ctx.rotate(73*5*Math.PI/180);
ctx.translate(-810,-620);


ctx.translate(810,649);
ctx.rotate(-b*Math.PI/180);
ctx.strokeStyle="#2E9AFE";
ctx.beginPath();
ctx.moveTo(120,0);
ctx.lineTo(130,-10);
ctx.lineTo(130,10);
ctx.lineTo(120,0);
//ctx.lineWidth=2;
ctx.strokeStyle="#2E9AFE";
ctx.closePath();
ctx.stroke();

ctx.rotate(b*Math.PI/180);
ctx.translate(-810,-649);

ctx.translate(660,470)
ctx.beginPath();
ctx.strokeStyle="yellow";
ctx.lineWidth=2;
ctx.moveTo(140,150);
ctx.lineTo(160,150);
ctx.moveTo(150,160);
ctx.lineTo(150,143);
ctx.moveTo(147,158);
ctx.lineTo(153,158);
ctx.closePath();
ctx.stroke();
ctx.beginPath();
ctx.moveTo(150,50);
ctx.lineTo(150,70)
ctx.lineWidth=2;
ctx.strokeStyle="yellow";
ctx.closePath();
ctx.stroke();


ctx.font="10px Arial";
ctx.fillStyle = "blue";
ctx.textAlign = "left";
ctx.fillText("10",105,185);
ctx.fillText("20",80,200);

ctx.beginPath();
ctx.strokeStyle="#80FF00";
ctx.moveTo(150,143);
ctx.lineTo(150,80);
ctx.lineTo(146,72.5);
ctx.lineTo(150,65);
ctx.lineTo(154,72.5);
ctx.lineTo(150,80);
ctx.closePath();
ctx.stroke();


ctx.font="8px Arial";
ctx.fillStyle = "purple";
ctx.textAlign = "left";
ctx.fillText("*EDDL",150,220);
ctx.translate(-660,-470);
              }
