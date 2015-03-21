/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

var gammln = function(xx) {
    var x,y,tmp,ser, j;
    var cof = [76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5];

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*Math.log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+Math.log(2.5066282746310005*ser/x);
}

var ITMAX = 100;
var EPS = 3.0e-7;
var FPMIN = 1.0e-30;

var gcf = function(a, x) {
    var i;
    var an,b,c,d,del,h;

    var gln=gammln(a);
    b=x+1.0-a;
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for (i=1;i<=ITMAX;i++) {
        an = -i*(i-a);
        b += 2.0;
        d=an*d+b;
        if (Math.abs(d) < FPMIN) d=FPMIN;
        c=b+an/c;
        if (Math.abs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (Math.abs(del-1.0) < EPS) break;
    }
    if (i > ITMAX) alert("a too large, ITMAX too small in gcf");

    return {
        gammcf: Math.exp(-x+a*Math.log(x)-(gln))*h,
        gln: gln
    };
}

var gser = function(a, x) {
    var n;
    var sum,del,ap;

    var gln=gammln(a), gamser;
    if (x <= 0.0) {
        if (x < 0.0) alert("x less than 0 in routine gser");
        gamser=0.0;
    } else {
        ap=a;
        del=sum=1.0/a;
        for (n=1;n<=ITMAX;n++) {
            ++ap;
            del *= x/ap;
            sum += del;
            if (Math.abs(del) < Math.abs(sum)*EPS) {
                gamser=sum*Math.exp(-x+a*Math.log(x)-(gln));
                return { gamser: gamser, gln: gln, };
            }
        }
        alert("a too large, ITMAX too small in routine gser");
    }
    return { gamser: gamser, gln: gln, };
}

var gammp = function(a, x) {
    if (x < 0.0 || a <= 0.0) alert("Invalid arguments in routine gammp");
    if (x < (a+1.0)) {
        var ret = gser(a,x);
        return ret.gamser;
    } else {
        var ret = gcf(a,x);
        return 1.0 - ret.gammcf;
    }
}

var erf = function(x) {
    // constants
    var a1 =  0.254829592;
    var a2 = -0.284496736;
    var a3 =  1.421413741;
    var a4 = -1.453152027;
    var a5 =  1.061405429;
    var p  =  0.3275911;

    // Save the sign of x
    var sign = 1;
    if (x < 0) {
        sign = -1;
    }
    x = Math.abs(x);

    // A&S formula 7.1.26
    var t = 1.0/(1.0 + p*x);
    var y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*Math.exp(-x*x);

    return sign*y;
}

