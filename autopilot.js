/**
 * @author mrdoob / http://mrdoob.com/
 * @author alteredq / http://alteredqualia.com/
 * @author paulirish / http://paulirish.com/
 */

THREE.FirstPersonControls = function ( object, domElement ) {

	this.object = object;
	this.target = new THREE.Vector3( 0, 0, 0 );

	this.domElement = ( domElement !== undefined ) ? domElement : document;

	this.enabled = true;

	this.movementSpeed = 1.0;
	this.lookSpeed = 0.005;

	this.lookVertical = true;
	this.autoForward = false;

	this.activeLook = true;

	this.heightSpeed = false;
	this.heightCoef = 1.0;
	this.heightMin = 0.0;
	this.heightMax = 1.0;

	this.constrainVertical = false;
	this.verticalMin = 0;
	this.verticalMax = Math.PI;

	this.autoSpeedFactor = 0.0;

	this.mouseX = 0;
	this.mouseY = 0;

	this.lat = 0;
	this.lon = 0;
	this.phi = 0;
	this.theta = 0;

	this.moveForward = false;
	this.moveBackward = false;
	this.moveLeft = false;
	this.moveRight = false;

	this.mouseDragOn = false;

	this.viewHalfX = 0;
	this.viewHalfY = 0;

	if ( this.domElement !== document ) {

		this.domElement.setAttribute( 'tabindex', - 1 );

	}

	//

	this.handleResize = function () {

		if ( this.domElement === document ) {

			this.viewHalfX = window.innerWidth / 2;
			this.viewHalfY = window.innerHeight / 2;

		} else {

			this.viewHalfX = this.domElement.offsetWidth / 2;
			this.viewHalfY = this.domElement.offsetHeight / 2;

		}

	};
/*
	this.onMouseDown = function ( event ) {

		if ( this.domElement !== document ) {

			this.domElement.focus();

		}

		event.preventDefault();
		event.stopPropagation();

		if ( this.activeLook ) {

			switch ( event.button ) {

				case 0: this.moveForward = true; break;
				case 2: this.moveBackward = true; break;

			}

		}

		this.mouseDragOn = true;

	};

	this.onMouseUp = function ( event ) {

		event.preventDefault();
		event.stopPropagation();

		if ( this.activeLook ) {

			switch ( event.button ) {

				case 0: this.moveForward = false; break;
				case 2: this.moveBackward = false; break;

			}

		}

		this.mouseDragOn = false;

	};

	this.onMouseMove = function ( event ) {

		if ( this.domElement === document ) {

			this.mouseX = event.pageX - this.viewHalfX;
			this.mouseY = event.pageY - this.viewHalfY;

		} else {

			this.mouseX = event.pageX - this.domElement.offsetLeft - this.viewHalfX;
			this.mouseY = event.pageY - this.domElement.offsetTop - this.viewHalfY;

		}

	};

	this.onKeyDown = function ( event ) {

		//event.preventDefault();

		switch ( event.keyCode ) {

			case 38: 
			case 87:  this.moveForward = true; break;

			case 37: 
			case 65:  this.moveLeft = true; break;

			case 40: 
			case 83:  this.moveBackward = true; break;

			case 39: 
			case 68:  this.moveRight = true; break;

			case 82:  this.moveUp = true; break;
			case 70:  this.moveDown = true; break;

		}

	};

	this.onKeyUp = function ( event ) {

		switch ( event.keyCode ) {

			case 38: 
			case 87:  this.moveForward = false; break;

			case 37: 
			case 65:  this.moveLeft = false; break;

			case 40: 
			case 83:  this.moveBackward = false; break;

			case 39: 
			case 68:  this.moveRight = false; break;

			case 82:  this.moveUp = false; break;
			case 70:  this.moveDown = false; break;

		}

	};*/
  aero = function () 
{

    var S = 0.017;
    var AR = 0.86;
    var e = 0.9;
    var m = 0.003;
    var g = 9.8;
    var rho = 1.225;
    var CLa = Math.PI * AR / (1 + Math.sqrt(1 + Math.pow((AR / 2), 2)));
    var CDo = 0.02;
    var epsilon = 1 / (Math.PI * e * AR);
    var CL = Math.sqrt(CDo / epsilon);
    var CD = CDo + epsilon * Math.pow(CL, 2);
    var LDmax = CL / CD;
    var Gam = -Math.atan(1 / LDmax);
    var V = Math.sqrt(2 * m * g / (rho * S * (CL * Math.cos(Gam) - CD * Math.sin(Gam))));
    var Alpha = CL / CLa;
    var q = 0.5 * rho * Math.pow(V, 2);

    var H = 2;
    var R = 0;
    var to = 0;
    var tf = 6;
    var xo = [];
    xo[0] = 3*V;
    xo[1] = Gam;
    xo[2] = H;
    xo[3] = R;

    var stepsize = 1 / 3600;
    var s1 = [];
    var s2 = [];
    var s3 = [];
    var tn = [];
    var xn = [];
    for (var i = 0; i < 4; i++)
    {
        xn[i] = [];
    }

    tn[0] = to;
    for (var i = 0; i < 4; i++)
    {
        xn[i][0] = xo[i];
    }

    for (var i = 1; i <= (1 / stepsize * (tf - to) + 1); i++)
    {
        tn[i] = to + i * stepsize;

        q = 0.5 * rho * Math.pow(xn[0][i - 1], 2);

        s1[0] = (-1 * CD * q * S - m * g * Math.sin(xn[1][i - 1])) / m;
        s1[1] = (CL * q * S - m * g * Math.cos(xn[1][i - 1])) / (m * xn[0][i - 1]);
        s1[2] = xn[0][i - 1] * Math.sin(xn[1][i - 1]);
        s1[3] = xn[0][i - 1] * Math.cos(xn[1][i - 1]);

        s2[0] = (-1 * CD * q * S - m * g * Math.sin(xn[1][i - 1] + (stepsize / 2) * s1[1])) / m;
        s2[1] = (CL * q * S - m * g * Math.cos(xn[1][i - 1] + (stepsize / 2) * s1[1])) / (m * (xn[0][i - 1] + (stepsize / 2) * s1[0]));
        s2[2] = (xn[0][i - 1] + (stepsize / 2) * s1[0]) * Math.sin(xn[1][i - 1] + (stepsize / 2) * s1[1]);
        s2[3] = (xn[0][i - 1] + (stepsize / 2) * s1[0]) * Math.cos(xn[1][i - 1] + (stepsize / 2) * s1[1]);

        s3[0] = (-1 * CD * q * S - m * g * Math.sin(xn[1][i - 1] + (stepsize / 4) * 3 * s2[1])) / m;
        s3[1] = (CL * q * S - m * g * Math.cos(xn[1][i - 1] + (stepsize / 4) * 3 * s2[1])) / (m * (xn[0][i - 1] + (stepsize / 4) * 3 * s2[0]));
        s3[2] = (xn[0][i - 1] + (stepsize / 4) * 3 * s2[0]) * Math.sin((xn[1][i - 1] + (stepsize / 4) * 3 * s2[1]));
        s3[3] = (xn[0][i - 1] + (stepsize / 4) * 3 * s2[0]) * Math.cos((xn[1][i - 1] + (stepsize / 4) * 3 * s2[1]));

        for (var j = 0; j < 4; ++j)
        {
            xn[j][i] = xn[j][i - 1] + (stepsize / 9) * (2 * s1[j] + 3 * s2[j] + 4 * s3[j]);
        }
    }

    var Vn = [];
    for (var i = 0; i <= 1 / stepsize * (tf - to) + 1; i++)
    {
        Vn[i] = xn[0][i];
    }
    var Gamn = [];
    for (var i = 0; i <= 1 / stepsize * (tf - to) + 1; i++)
    {
        Gamn[i] = xn[1][i];
    }
    var Hn = [];
    for (var i = 0; i <= 1 / stepsize * (tf - to) + 1; i++)
    {
        Hn[i] = xn[2][i];
    }
    var Rn = [];
    for (var i = 0; i <= 1 / stepsize * (tf - to) + 1; i++)
    {
        Rn[i] = xn[3][i];
    }
};

	this.update = function( delta ) {

		if ( this.enabled === false ) return;

		if ( this.heightSpeed ) {

			var y = THREE.Math.clamp( this.object.position.y, this.heightMin, this.heightMax );
			var heightDelta = y - this.heightMin;

			this.autoSpeedFactor = delta * ( heightDelta * this.heightCoef );

		} else {

			this.autoSpeedFactor = 0.0;

		}

		var actualMoveSpeed = delta * this.movementSpeed;

		if ( this.moveForward || ( this.autoForward && ! this.moveBackward ) ) this.object.translateZ( - ( actualMoveSpeed + this.autoSpeedFactor ) );
		if ( this.moveBackward ) this.object.translateZ( actualMoveSpeed );

		if ( this.moveLeft ) this.object.translateX( - actualMoveSpeed );
		if ( this.moveRight ) this.object.translateX( actualMoveSpeed );

		if ( this.moveUp ) this.object.translateY( actualMoveSpeed );
		if ( this.moveDown ) this.object.translateY( - actualMoveSpeed );

		var actualLookSpeed = delta * this.lookSpeed;

		if ( ! this.activeLook ) {

			actualLookSpeed = 0;

		}

		var verticalLookRatio = 1;

		if ( this.constrainVertical ) {

			verticalLookRatio = Math.PI / ( this.verticalMax - this.verticalMin );

		}

		this.lon += this.mouseX * actualLookSpeed;
		if ( this.lookVertical ) this.lat -= this.mouseY * actualLookSpeed * verticalLookRatio;

		this.lat = Math.max( - 85, Math.min( 85, this.lat ) );
		this.phi = THREE.Math.degToRad( 90 - this.lat );

		this.theta = THREE.Math.degToRad( this.lon );

		if ( this.constrainVertical ) {

			this.phi = THREE.Math.mapLinear( this.phi, 0, Math.PI, this.verticalMin, this.verticalMax );

		}

		var targetPosition = this.target,
			position = this.object.position;

		targetPosition.x = position.x + 100 * Math.sin( this.phi ) * Math.cos( this.theta );
		targetPosition.y = position.y + 100 * Math.cos( this.phi );
		targetPosition.z = position.z + 100 * Math.sin( this.phi ) * Math.sin( this.theta );

		this.object.lookAt( targetPosition );

	};

	function contextmenu( event ) {

		event.preventDefault();

	}

	this.dispose = function() {

		this.domElement.removeEventListener( 'contextmenu', contextmenu, false );
		this.domElement.removeEventListener( 'mousedown', _onMouseDown, false );
		this.domElement.removeEventListener( 'mousemove', _onMouseMove, false );
		this.domElement.removeEventListener( 'mouseup', _onMouseUp, false );

		window.removeEventListener( 'keydown', _onKeyDown, false );
		window.removeEventListener( 'keyup', _onKeyUp, false );

	};

	var _onMouseMove = bind( this, this.onMouseMove );
	var _onMouseDown = bind( this, this.onMouseDown );
	var _onMouseUp = bind( this, this.onMouseUp );
	var _onKeyDown = bind( this, this.onKeyDown );
	var _onKeyUp = bind( this, this.onKeyUp );

	this.domElement.addEventListener( 'contextmenu', contextmenu, false );
	this.domElement.addEventListener( 'mousemove', _onMouseMove, false );
	this.domElement.addEventListener( 'mousedown', _onMouseDown, false );
	this.domElement.addEventListener( 'mouseup', _onMouseUp, false );

	window.addEventListener( 'keydown', _onKeyDown, false );
	window.addEventListener( 'keyup', _onKeyUp, false );

	function bind( scope, fn ) {

		return function () {

			fn.apply( scope, arguments );

		};

	}

	this.handleResize();

};
