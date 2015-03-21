var vec = function(len, def) {
    var vec = [];
    for(var i = 0; i < len; i++) {
        vec.push(def);
    }
    return vec;
}
var mat = function(x, y, def) {
    if(def == null) {
        def = 0;
    }
    var mat = [];
    for(var i = 0; i < x; i++) {
        mat.push(vec(y, def));
    }
    return mat;
}

var fx = function(xx, sx3, xcenter, sx) {
    return sx3 * Math.tan((xx - xcenter) * sx / sx3);
}

var debye = function(scale, rg, bkg, x) {
    var res = x.concat([]);
    for(var i = 0; i < res.length; i++) {
        var v = res[i];
        var qr2 = (v * rg) * (v * rg);
        var pq = 2 * (Math.exp(-qr2) - 1 + qr2)/(qr2*qr2);
        pq *= scale;
        res[i] = pq + bkg;
    }
    return res;
}

var increment_pixel = function(pixel, ddr, dxx, dyy, aveint, dsq, ncells, nq, nd2) {
    var ir = (0|(Math.sqrt(dxx*dxx + dyy*dyy)/ddr)) + 1;
    if(ir > nq) {
        nq = ir; // ORIG: reset maximum number of q-values
    }
    aveint[ir-1] += pixel / nd2;
    dsq[ir-1] += pixel * pixel / nd2;
    ncells[ir-1] += 1 / nd2;

    return nq;
}

var get_resolution = function(inq, lambda, lambdaWidth, ddet, apoff,
        s1, s2, l1, l2, bs, del_r) {
    var vz_1 = 3.956e5; // ORIG: velocity [cm/s] of 1 A neutron
    var g = 981.0; // ORIG: gravity acceleration [cm/s^2]

    var use_lense = false; //TODO:

    s1 *= 0.5 * 0.1; // ORIG: convert to radisu and [cm]
    s2 *= 0.5 * 0.1;
    l1 *= 100; // ORIG: [cm]
    l1 -= apoff; // ORIG: correct the distance

    l2 *= 100;
    l2 += apoff;

    bs *= 0.5 * 0.1; // ORIG: convert to radisu and [cm]
    del_r *= 0.1; // ORIG: width of annulus, convert mm to [cm]

    // ORIG: start resolution calculation
    var a2 = s1*l2/l1 + s2*(l1+l2)/l1;
    var q_small = 2.0 * Math.PI * (bs - a2) * (1.0 - lambdaWidth) / (lambda * l2);
    var lp = 1.0 / (1.0/l1 + 1.0/l2);
    var v_lambda = lambdaWidth * lambdaWidth / 6.0;

    var v_b;
    if(use_lense) {
        v_b = 0.25 * (s1 * l2/l1) * (s1 * l2/l1) +
                0.25 * 2/3 * (lambdaWidth/lambda) * (lambdaWidth/lambda) *
                (s2*l2/lp) * (s2*l2/lp);
    } else {
        v_b = 0.25 * (s1 * l2/l1) * (s1 * l2/l1) +
                0.25 * (s2*l2/lp) * (s2*l2/lp);
    }
    var v_d = (ddet / 2.3548) * (ddet / 2.3548) + del_r * del_r / 12.0;
    var vz = vz_1 / lambda;
    var yg = 0.5 * g * l2 * (l1+l2)/(vz*vz);
    var v_g = 2.0*(2.0 * yg * yg * v_lambda);
    var r0 = l2 * Math.tan(2.0 * Math.asin(lambda * inq / (4.0 * Math.PI)));
    var delta = 0.5 * (bs - r0) * (bs - r0) / v_d;

    var inc_gamma;
    if(r0 < bs) {
        inc_gamma = Math.exp(gammln(1.5)) * (1 - gammp(1.5, delta));
    } else {
        inc_gamma = Math.exp(gammln(1.5)) * (1 + gammp(1.5, delta));
    }

    var fsubs = 0.5 * (1.0 + erf( (r0-bs) / Math.sqrt(2.0*v_d)));
    if(fsubs <= 0) {
        fsubs = 1e-10;
    }
    var fr = 1.0 + Math.sqrt(v_d) * Math.exp(-1.0 * delta) /
            (r0 * fsubs * Math.sqrt(2.0 * Math.PI));
    var fv = inc_gamma / (fsubs * Math.sqrt(Math.PI)) - r0*r0 * (fr-1.0)*(fr-1.0)/v_d;

    var rmd = fr * r0;
    var v_r1 = v_b + fv*v_d + v_g;
    var rm = rmd + 0.5*v_r1/rmd;
    var v_r = v_r1 - 0.5*(v_r1/rmd)*(v_r1/rmd);
    if(v_r < 0.0) {
        v_r = 0.0;
    }

    var qbar = (4.0 * Math.PI / lambda) * Math.sin(0.5 * Math.atan(rm/l2));
    var sigmaq = qbar * Math.sqrt(v_r/(rmd*rmd) + v_lambda);

    return {
        sigmaq: sigmaq,
        qbar: qbar,
        fsubs: fsubs,
    };
}

// 1, 54.8, 5
var sourceToSampleDist = function(ng, s12, l2diff, type) {
    if(type == "huber") {
        return (1632 - 155 * ng - s12 - l2diff) / 100.0;
    } else if(type == "chamber") {
        return (1632 - 155 * ng - l2diff) / 100.0;
    }
    alert("ERROR: sourceToSampleDist");
}

var Stat = function() {
    var stat = this;
    stat.pixelsX = 128;
    stat.pixelsY = 128;
    stat.xcenter = stat.pixelsX/2 + 0.5;
    stat.ycenter = stat.pixelsY/2 + 0.5;

    stat.x0 = 64;
    stat.y0 = 64;
    stat.sx = 5.08; // mm/pixel(x)
    stat.sx3 = 10000; //nonlinear coeff
    stat.sy = 5.08; // mm/pixel(y)
    stat.sy3 = 10000; //nonlinear coeff

    stat.dtsize = 650; // detector size in mm
    stat.dtdist = 13170; // detector distance in mm

    stat.dr = 1; // annulus width set by user, default is one
    stat.ddr = stat.dr * stat.sx; //step size, in mm

    stat.rcentr = 100;
    stat.large_num = 1;
    stat.small_num = 1e-10;

    stat.bs = 76.2;
    stat.s1 = 50;
    stat.s2 = 6.35;
    stat.l1 = sourceToSampleDist(1, 54.8, 5, "chamber");
    stat.l2 = 13.17;
    stat.lambdaWidth = 0.125;
    stat.apoff = 5.0;
    stat.ddet = 0.5;

    stat.lambda = 8.4;
};

Stat.prototype.calculate = function() {
    var stat = this;

    var qval = vec(500, 0);
    var aveint = vec(500, 0);
    var ncells = vec(500, 0);
    var dsq = vec(500, 0);
    var sigave = vec(500, 0);
    var qbar = vec(500, 0);
    var sigmaq = vec(500, 0);
    var fsubs = vec(500, 0);

    var dxbm = fx(stat.x0, stat.sx3, stat.xcenter, stat.sx);
    var dybm = fx(stat.y0, stat.sy3, stat.ycenter, stat.sy);

    var nq = 1;

    var data = mat(stat.pixelsX, stat.pixelsY, 1);

    // mask three pixels all around, I don't know why...
    // many loops are based on 1-based index, because original subroutine was
    // written in FORTRAN, sigh...
    for(var i = 3; i < stat.pixelsX-3; i++) {
        var dxi = fx(i+1, stat.sx3, stat.xcenter, stat.sx);
        var dx = dxi - dxbm;
        for(var j = 3; j < stat.pixelsY-3; j++) {
            var dyi = fx(j+1, stat.sy3, stat.ycenter, stat.sy);
            var dy = dyi - dybm;

            var dr2 = Math.sqrt(dx*dx + dy*dy);

            var nd, fd;
            if(dr2 > stat.rcentr) {
                // ORIG: keep pixel whole
                nd = 1;
                fd = 1;
            } else {
                // ORIG: break pixel into 9 equal parts
                nd = 3;
                fd = 2;
            }
            var nd2 = nd * nd;
            for(var l = 1; l <= nd; l++) {
                var dxx = dx + (l - fd)*stat.sx/3;
                for(var k = 1; k <= nd; k++) {
                    var dyy = dy + (k - fd)*stat.sy/3;

                    nq = increment_pixel(data[i][j], stat.ddr, dxx, dyy,
                        aveint, dsq, ncells, nq, nd2);
                }
            }
        }
    }

    var ntotal = 0;
    for(var i = 0; i < nq; i++) {
        var rr = (2*i + 1) * stat.ddr / 2;
        var theta = 0.5 * Math.atan(rr / stat.dtdist);
        qval[i] = (4 * Math.PI / stat.lambda) * Math.sin(theta);
        if(ncells[i] == 0) {
            // ORIG: no pixels in annuli, data unknown
            aveint[i] = 0;
            sigave[i] = stat.large_num;
        } else {
            if(ncells[i] <= 1) {
                // ORIG: need more than one pixel to determine error
                aveint[i] /= ncells[i];
                sigave[i] = stat.large_num;
            } else {
                // ORIG: assume that the intensity in each pixel in annuli is normally
                // ORIG: distributed about mean...
                aveint[i] /= ncells[i];
                var avesq = (aveint[i] * aveint[i]);
                var aveisq = dsq[i] / ncells[i];
                var v = aveisq - avesq;
                if(v <= 0) {
                    sigave[i] = stat.small_num;
                } else {
                    sigave[i] = Math.sqrt(v / (ncells[i] - 1));
                }
            }
        }
        ntotal += ncells[i];
    }

    // ORIG: Do the extra 3 columns of resolution calculations starting here.
    for(var i = 0; i < nq; i++) {
        var obj = get_resolution(qval[i],
            stat.lambda, stat.lambdaWidth,
            stat.ddet, stat.apoff,
            stat.s1, stat.s2, stat.l1, stat.l2, stat.bs, stat.ddr);

        sigmaq[i] = obj.sigmaq;
        qbar[i] = obj.qbar;
        fsubs[i] = obj.fsubs;
    }
    for(var i = 0; i < fsubs.length; i++) {
        // ORIG: keep the values from being too small
        fsubs[i] += 1e-8;
    }

    aveint = debye(1000, 100, 0.0, qval);
    for(var i = 0; i < aveint.length; i++) {
        aveint[i] *= fsubs[i];
    }

    var len = qval.indexOf(0);
    qval.splice(len);
    aveint.splice(len);
    return {
        qval: qval,
        aveint: aveint,
    };
};

