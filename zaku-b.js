 
window.zaku = function(){


  this.utils = {

     keyframe : [0,1],
     newValue : null,
     initialValue : null

  }

};



var pathSegmentPattern = /[a-z][^a-z]*/ig;

//Starting of Arc to Bezier Utils
var TAU = Math.PI * 2;


/* eslint-disable space-infix-ops */

// Calculate an angle between two unit vectors
//
// Since we measure angle between radii of circular arcs,
// we can use simplified math (without length normalization)
//
function unit_vector_angle(ux, uy, vx, vy) {
  var sign = (ux * vy - uy * vx < 0) ? -1 : 1;
  var dot  = ux * vx + uy * vy;

  // Add this to work with arbitrary vectors:
  // dot /= Math.sqrt(ux * ux + uy * uy) * Math.sqrt(vx * vx + vy * vy);

  // rounding errors, e.g. -1.0000000000000002 can screw up this
  if (dot >  1.0) { dot =  1.0; }
  if (dot < -1.0) { dot = -1.0; }

  return sign * Math.acos(dot);
}


// Convert from endpoint to center parameterization,
// see http://www.w3.org/TR/SVG11/implnote.html#ArcImplementationNotes
//
// Return [cx, cy, theta1, delta_theta]
//
function get_arc_center(x1, y1, x2, y2, fa, fs, rx, ry, sin_phi, cos_phi) {
  // Step 1.
  //
  // Moving an ellipse so origin will be the middlepoint between our two
  // points. After that, rotate it to line up ellipse axes with coordinate
  // axes.
  //
  var x1p =  cos_phi*(x1-x2)/2 + sin_phi*(y1-y2)/2;
  var y1p = -sin_phi*(x1-x2)/2 + cos_phi*(y1-y2)/2;

  var rx_sq  =  rx * rx;
  var ry_sq  =  ry * ry;
  var x1p_sq = x1p * x1p;
  var y1p_sq = y1p * y1p;

  // Step 2.
  //
  // Compute coordinates of the centre of this ellipse (cx', cy')
  // in the new coordinate system.
  //
  var radicant = (rx_sq * ry_sq) - (rx_sq * y1p_sq) - (ry_sq * x1p_sq);

  if (radicant < 0) {
    // due to rounding errors it might be e.g. -1.3877787807814457e-17
    radicant = 0;
  }

  radicant /=   (rx_sq * y1p_sq) + (ry_sq * x1p_sq);
  radicant = Math.sqrt(radicant) * (fa === fs ? -1 : 1);

  var cxp = radicant *  rx/ry * y1p;
  var cyp = radicant * -ry/rx * x1p;

  // Step 3.
  //
  // Transform back to get centre coordinates (cx, cy) in the original
  // coordinate system.
  //
  var cx = cos_phi*cxp - sin_phi*cyp + (x1+x2)/2;
  var cy = sin_phi*cxp + cos_phi*cyp + (y1+y2)/2;

  // Step 4.
  //
  // Compute angles (theta1, delta_theta).
  //
  var v1x =  (x1p - cxp) / rx;
  var v1y =  (y1p - cyp) / ry;
  var v2x = (-x1p - cxp) / rx;
  var v2y = (-y1p - cyp) / ry;

  var theta1 = unit_vector_angle(1, 0, v1x, v1y);
  var delta_theta = unit_vector_angle(v1x, v1y, v2x, v2y);

  if (fs === 0 && delta_theta > 0) {
    delta_theta -= TAU;
  }
  if (fs === 1 && delta_theta < 0) {
    delta_theta += TAU;
  }

  return [ cx, cy, theta1, delta_theta ];
}

//
// Approximate one unit arc segment with bézier curves,
// see http://math.stackexchange.com/questions/873224
//
function approximate_unit_arc(theta1, delta_theta) {
  var alpha = 4/3 * Math.tan(delta_theta/4);

  var x1 = Math.cos(theta1);
  var y1 = Math.sin(theta1);
  var x2 = Math.cos(theta1 + delta_theta);
  var y2 = Math.sin(theta1 + delta_theta);

  return [ x1, y1, x1 - y1*alpha, y1 + x1*alpha, x2 + y2*alpha, y2 - x2*alpha, x2, y2 ];
}

function a2c(x1, y1, x2, y2, fa, fs, rx, ry, phi) {
  var sin_phi = Math.sin(phi * TAU / 360);
  var cos_phi = Math.cos(phi * TAU / 360);

  // Make sure radii are valid
  //
  var x1p =  cos_phi*(x1-x2)/2 + sin_phi*(y1-y2)/2;
  var y1p = -sin_phi*(x1-x2)/2 + cos_phi*(y1-y2)/2;

  if (x1p === 0 && y1p === 0) {
    // we're asked to draw line to itself
    return [];
  }

  if (rx === 0 || ry === 0) {
    // one of the radii is zero
    return [];
  }


  // Compensate out-of-range radii
  //
  rx = Math.abs(rx);
  ry = Math.abs(ry);

  var lambda = (x1p * x1p) / (rx * rx) + (y1p * y1p) / (ry * ry);
  if (lambda > 1) {
    rx *= Math.sqrt(lambda);
    ry *= Math.sqrt(lambda);
  }


  // Get center parameters (cx, cy, theta1, delta_theta)
  //
  var cc = get_arc_center(x1, y1, x2, y2, fa, fs, rx, ry, sin_phi, cos_phi);

  var result = [];
  var theta1 = cc[2];
  var delta_theta = cc[3];

  // Split an arc to multiple segments, so each segment
  // will be less than τ/4 (= 90°)
  //
  var segments = Math.max(Math.ceil(Math.abs(delta_theta) / (TAU / 4)), 1);
  delta_theta /= segments;

  for (var i = 0; i < segments; i++) {
    result.push(approximate_unit_arc(theta1, delta_theta));
    theta1 += delta_theta;
  }

  // We have a bezier approximation of a unit circle,
  // now need to transform back to the original ellipse
  //
  return result.map(function (curve) {
    for (var i = 0; i < curve.length; i += 2) {
      var x = curve[i + 0];
      var y = curve[i + 1];

      // scale
      x *= rx;
      y *= ry;

      // rotate
      var xp = cos_phi*x - sin_phi*y;
      var yp = sin_phi*x + cos_phi*y;

      // translate
      curve[i + 0] = xp + cc[0];
      curve[i + 1] = yp + cc[1];
    }

    return curve;
  });
};

//starting of Zakarc parser

var paramCounts = { a: 7, c: 6, h: 1, l: 2, m: 2, r: 4, q: 4, s: 4, t: 2, v: 1, z: 0 };

var SPECIAL_SPACES = [
  0x1680, 0x180E, 0x2000, 0x2001, 0x2002, 0x2003, 0x2004, 0x2005, 0x2006,
  0x2007, 0x2008, 0x2009, 0x200A, 0x202F, 0x205F, 0x3000, 0xFEFF
];

function isSpace(ch) {
  return (ch === 0x0A) || (ch === 0x0D) || (ch === 0x2028) || (ch === 0x2029) || // Line terminators
    // White spaces
    (ch === 0x20) || (ch === 0x09) || (ch === 0x0B) || (ch === 0x0C) || (ch === 0xA0) ||
    (ch >= 0x1680 && SPECIAL_SPACES.indexOf(ch) >= 0);
}

function isCommand(code) {
  /*eslint-disable no-bitwise*/
  switch (code | 0x20) {
    case 0x6D/* m */:
    case 0x7A/* z */:
    case 0x6C/* l */:
    case 0x68/* h */:
    case 0x76/* v */:
    case 0x63/* c */:
    case 0x73/* s */:
    case 0x71/* q */:
    case 0x74/* t */:
    case 0x61/* a */:
    case 0x72/* r */:
      return true;
  }
  return false;
}

function isDigit(code) {
  return (code >= 48 && code <= 57);   // 0..9
}

function isDigitStart(code) {
  return (code >= 48 && code <= 57) || /* 0..9 */
          code === 0x2B || /* + */
          code === 0x2D || /* - */
          code === 0x2E;   /* . */
}


function State(path) {
  this.index  = 0;
  this.path   = path;
  this.max    = path.length;
  this.result = [];
  this.param  = 0.0;
  this.err    = '';
  this.segmentStart = 0;
  this.data   = [];
}

function skipSpaces(state) {
  while (state.index < state.max && isSpace(state.path.charCodeAt(state.index))) {
    state.index++;
  }
}


function scanParam(state) {
  var start = state.index,
      index = start,
      max = state.max,
      zeroFirst = false,
      hasCeiling = false,
      hasDecimal = false,
      hasDot = false,
      ch;

  if (index >= max) {
    state.err = 'Zakarc: missed param (at pos ' + index + ')';
    return;
  }
  ch = state.path.charCodeAt(index);

  if (ch === 0x2B/* + */ || ch === 0x2D/* - */) {
    index++;
    ch = (index < max) ? state.path.charCodeAt(index) : 0;
  }

  // This logic is shamelessly borrowed from Esprima
  // https://github.com/ariya/esprimas
  //
  if (!isDigit(ch) && ch !== 0x2E/* . */) {
    state.err = 'Zaku: param should start with 0..9 or `.` (at pos ' + index + ')';
    return;
  }

  if (ch !== 0x2E/* . */) {
    zeroFirst = (ch === 0x30/* 0 */);
    index++;

    ch = (index < max) ? state.path.charCodeAt(index) : 0;

    if (zeroFirst && index < max) {
      // decimal number starts with '0' such as '09' is illegal.
      if (ch && isDigit(ch)) {
        state.err = 'Zaku: numbers started with `0` such as `09` are ilegal (at pos ' + start + ')';
        return;
      }
    }

    while (index < max && isDigit(state.path.charCodeAt(index))) {
      index++;
      hasCeiling = true;
    }
    ch = (index < max) ? state.path.charCodeAt(index) : 0;
  }

  if (ch === 0x2E/* . */) {
    hasDot = true;
    index++;
    while (isDigit(state.path.charCodeAt(index))) {
      index++;
      hasDecimal = true;
    }
    ch = (index < max) ? state.path.charCodeAt(index) : 0;
  }

  if (ch === 0x65/* e */ || ch === 0x45/* E */) {
    if (hasDot && !hasCeiling && !hasDecimal) {
      state.err = 'Zaku: invalid float exponent (at pos ' + index + ')';
      return;
    }

    index++;

    ch = (index < max) ? state.path.charCodeAt(index) : 0;
    if (ch === 0x2B/* + */ || ch === 0x2D/* - */) {
      index++;
    }
    if (index < max && isDigit(state.path.charCodeAt(index))) {
      while (index < max && isDigit(state.path.charCodeAt(index))) {
        index++;
      }
    } else {
      state.err = 'Zaku: invalid float exponent (at pos ' + index + ')';
      return;
    }
  }

  state.index = index;
  state.param = parseFloat(state.path.slice(start, index)) + 0.0;
}


function finalizeSegment(state) {
  var cmd, cmdLC;

  // Process duplicated commands (without comand name)

  // This logic is shamelessly borrowed from Raphael
  // https://github.com/DmitryBaranovskiy/raphael/
  //
  cmd   = state.path[state.segmentStart];
  cmdLC = cmd.toLowerCase();

  var params = state.data;

  if (cmdLC === 'm' && params.length > 2) {
    state.result.push([ cmd, params[0], params[1] ]);
    params = params.slice(2);
    cmdLC = 'l';
    cmd = (cmd === 'm') ? 'l' : 'L';
  }

  if (cmdLC === 'r') {
    state.result.push([ cmd ].concat(params));
  } else {

    while (params.length >= paramCounts[cmdLC]) {
      state.result.push([ cmd ].concat(params.splice(0, paramCounts[cmdLC])));
      if (!paramCounts[cmdLC]) {
        break;
      }
    }
  }
}


function scanSegment(state) {
  var max = state.max,
      cmdCode, comma_found, need_params, i;

  state.segmentStart = state.index;
  cmdCode = state.path.charCodeAt(state.index);

  if (!isCommand(cmdCode)) {
    state.err = 'Zaku: bad command ' + state.path[state.index] + ' (at pos ' + state.index + ')';
    return;
  }

  need_params = paramCounts[state.path[state.index].toLowerCase()];

  state.index++;
  skipSpaces(state);

  state.data = [];

  if (!need_params) {
    // Z
    finalizeSegment(state);
    return;
  }

  comma_found = false;

  for (;;) {
    for (i = need_params; i > 0; i--) {
      scanParam(state);
      if (state.err.length) {
        return;
      }
      state.data.push(state.param);

      skipSpaces(state);
      comma_found = false;

      if (state.index < max && state.path.charCodeAt(state.index) === 0x2C/* , */) {
        state.index++;
        skipSpaces(state);
        comma_found = true;
      }
    }

    // after ',' param is mandatory
    if (comma_found) {
      continue;
    }

    if (state.index >= state.max) {
      break;
    }

    // Stop on next segment
    if (!isDigitStart(state.path.charCodeAt(state.index))) {
      break;
    }
  }

  finalizeSegment(state);
}


/* Returns array of segments:
 *
 * [
 *   [ command, coord1, coord2, ... ]
 * ]
 */
function pathParse(svgPath) {
  var state = new State(svgPath);
  var max = state.max;

  skipSpaces(state);

  while (state.index < max && !state.err.length) {
    scanSegment(state);
  }

  if (state.err.length) {
    state.result = [];

  } else if (state.result.length) {

    if ('mM'.indexOf(state.result[0][0]) < 0) {
      state.err = 'Zaku: string should start with `M` or `m`';
      state.result = [];
    } else {
      state.result[0][0] = 'M';
    }
  }

  return {
    err: state.err,
    segments: state.result
  };
};


// //////////////////////////END OF PATHPARSE//////////////////////


zaku.prototype = {
	
	get:function(id)
	{

       return document.getElementById(id);

	},
    path:function(elem)
	{

      var path = elem.getAttribute('d');

      // return path.replace(/,/g , ' ');
      return path;

	},
	normalize:function(pathArray)
	{

       for (var i = 0; i < pathArray.length; i++) 
        {
       	    var lineString = pathArray[i];

       	    var secondChar = lineString[1];
       	    if(secondChar === ' ')
       	    	pathArray[i] = pathArray[i].replace(' ' , '');
        }

       return pathArray;

	},

  zakarc:function(path) {
    
    var pstate = pathParse(path);

    // Array of path segments.
    // Each segment is array [command, param1, param2, ...]
    this.segments = pstate.segments;

    // Error message on parse error.
    this.err      = pstate.err;

    // Transforms stack for lazy evaluation
    this.__stack    = [];

    this.zakarc.prototype = {

         
       __matrix : function (m) {
          var self = this, i;

          // Quick leave for empty matrix
          if (!m.queue.length) { return; }

          this.iterate(function (s, index, x, y) {
            var p, result, name, isRelative;

            switch (s[0]) {

              // Process 'assymetric' commands separately
              case 'v':
                p      = m.calc(0, s[1], true);
                result = (p[0] === 0) ? [ 'v', p[1] ] : [ 'l', p[0], p[1] ];
                break;

              case 'V':
                p      = m.calc(x, s[1], false);
                result = (p[0] === m.calc(x, y, false)[0]) ? [ 'V', p[1] ] : [ 'L', p[0], p[1] ];
                break;

              case 'h':
                p      = m.calc(s[1], 0, true);
                result = (p[1] === 0) ? [ 'h', p[0] ] : [ 'l', p[0], p[1] ];
                break;

              case 'H':
                p      = m.calc(s[1], y, false);
                result = (p[1] === m.calc(x, y, false)[1]) ? [ 'H', p[0] ] : [ 'L', p[0], p[1] ];
                break;

              case 'a':
              case 'A':
                // ARC is: ['A', rx, ry, x-axis-rotation, large-arc-flag, sweep-flag, x, y]

                // Drop segment if arc is empty (end point === start point)
                /*if ((s[0] === 'A' && s[6] === x && s[7] === y) ||
                    (s[0] === 'a' && s[6] === 0 && s[7] === 0)) {
                  return [];
                }*/

                // Transform rx, ry and the x-axis-rotation
                var ma = m.toArray();
                var e = ellipse(s[1], s[2], s[3]).transform(ma);

                // flip sweep-flag if matrix is not orientation-preserving
                if (ma[0] * ma[3] - ma[1] * ma[2] < 0) {
                  s[5] = s[5] ? '0' : '1';
                }

                // Transform end point as usual (without translation for relative notation)
                p = m.calc(s[6], s[7], s[0] === 'a');

                // Empty arcs can be ignored by renderer, but should not be dropped
                // to avoid collisions with `S A S` and so on. Replace with empty line.
                if ((s[0] === 'A' && s[6] === x && s[7] === y) ||
                    (s[0] === 'a' && s[6] === 0 && s[7] === 0)) {
                  result = [ s[0] === 'a' ? 'l' : 'L', p[0], p[1] ];
                  break;
                }

                // if the resulting ellipse is (almost) a segment ...
                if (e.isDegenerate()) {
                  // replace the arc by a line
                  result = [ s[0] === 'a' ? 'l' : 'L', p[0], p[1] ];
                } else {
                  // if it is a real ellipse
                  // s[0], s[4] and s[5] are not modified
                  result = [ s[0], e.rx, e.ry, e.ax, s[4], s[5], p[0], p[1] ];
                }

                break;

              case 'm':
                // Edge case. The very first `m` should be processed as absolute, if happens.
                // Make sense for coord shift transforms.
                isRelative = index > 0;

                p = m.calc(s[1], s[2], isRelative);
                result = [ 'm', p[0], p[1] ];
                break;

              default:
                name       = s[0];
                result     = [ name ];
                isRelative = (name.toLowerCase() === name);

                // Apply transformations to the segment
                for (i = 1; i < s.length; i += 2) {
                  p = m.calc(s[i], s[i + 1], isRelative);
                  result.push(p[0], p[1]);
                }
            }

            zaku.segments[index] = result;
          }, true);
        },


      __evaluateStack : function () {
        var m, i;

        if (!zaku.__stack.length) { return; }

        if (zaku.__stack.length === 1) {
          zaku.__matrix(zaku.__stack[0]);
          zaku.__stack = [];
          return;
        }

        m = matrix();
        i = zaku.__stack.length;

        while (--i >= 0) {
          m.matrix(zaku.__stack[i].toArray());
        }

        zaku.__matrix(m);
        zaku.__stack = [];
      },

      matrix : function (m) {
        zaku.__stack.push(matrix().matrix(m));
        return zaku;
      },

     round : function (d) {
      var contourStartDeltaX = 0, contourStartDeltaY = 0, deltaX = 0, deltaY = 0, l;

      d = d || 0;

      this.__evaluateStack();

       zaku.segments.forEach(function (s) {
        var isRelative = (s[0].toLowerCase() === s[0]);

        switch (s[0]) {
          case 'H':
          case 'h':
            if (isRelative) { s[1] += deltaX; }
            deltaX = s[1] - s[1].toFixed(d);
            s[1] = +s[1].toFixed(d);
            return;

          case 'V':
          case 'v':
            if (isRelative) { s[1] += deltaY; }
            deltaY = s[1] - s[1].toFixed(d);
            s[1] = +s[1].toFixed(d);
            return;

          case 'Z':
          case 'z':
            deltaX = contourStartDeltaX;
            deltaY = contourStartDeltaY;
            return;

          case 'M':
          case 'm':
            if (isRelative) {
              s[1] += deltaX;
              s[2] += deltaY;
            }

            deltaX = s[1] - s[1].toFixed(d);
            deltaY = s[2] - s[2].toFixed(d);

            contourStartDeltaX = deltaX;
            contourStartDeltaY = deltaY;

            s[1] = +s[1].toFixed(d);
            s[2] = +s[2].toFixed(d);
            return;

          case 'A':
          case 'a':
            // [cmd, rx, ry, x-axis-rotation, large-arc-flag, sweep-flag, x, y]
            if (isRelative) {
              s[6] += deltaX;
              s[7] += deltaY;
            }

            deltaX = s[6] - s[6].toFixed(d);
            deltaY = s[7] - s[7].toFixed(d);

            s[1] = +s[1].toFixed(d);
            s[2] = +s[2].toFixed(d);
            s[3] = +s[3].toFixed(d + 2); // better precision for rotation
            s[6] = +s[6].toFixed(d);
            s[7] = +s[7].toFixed(d);
            return;

          default:
            // a c l q s t
            l = s.length;

            if (isRelative) {
              s[l - 2] += deltaX;
              s[l - 1] += deltaY;
            }

            deltaX = s[l - 2] - s[l - 2].toFixed(d);
            deltaY = s[l - 1] - s[l - 1].toFixed(d);

            s.forEach(function (val, i) {
              if (!i) { return; }
              s[i] = +s[i].toFixed(d);
            });
            return;
        }
      });

      return zaku;
    },

   iterate : function (iterator, keepLazyStack) {
      var segments = zaku.segments,
          replacements = {},
          needReplace = false,
          lastX = 0,
          lastY = 0,
          countourStartX = 0,
          countourStartY = 0;
      var i, j, newSegments;

      if (!keepLazyStack) {
        this.__evaluateStack();
      }

      segments.forEach(function (s, index) {

        var res = iterator(s, index, lastX, lastY);

        if (Array.isArray(res)) {
          replacements[index] = res;
          needReplace = true;
        }

        var isRelative = (s[0] === s[0].toLowerCase());

        // calculate absolute X and Y
        switch (s[0]) {
          case 'm':
          case 'M':
            lastX = s[1] + (isRelative ? lastX : 0);
            lastY = s[2] + (isRelative ? lastY : 0);
            countourStartX = lastX;
            countourStartY = lastY;
            return;

          case 'h':
          case 'H':
            lastX = s[1] + (isRelative ? lastX : 0);
            return;

          case 'v':
          case 'V':
            lastY = s[1] + (isRelative ? lastY : 0);
            return;

          case 'z':
          case 'Z':
            // That make sence for multiple contours
            lastX = countourStartX;
            lastY = countourStartY;
            return;

          default:
            lastX = s[s.length - 2] + (isRelative ? lastX : 0);
            lastY = s[s.length - 1] + (isRelative ? lastY : 0);
        }
      });

      // Replace segments if iterator return results

      if (!needReplace) { return this; }

      newSegments = [];

      for (i = 0; i < segments.length; i++) {
        if (typeof replacements[i] !== 'undefined') {
          for (j = 0; j < replacements[i].length; j++) {
            newSegments.push(replacements[i][j]);
          }
        } else {
          newSegments.push(segments[i]);
        }
      }

      zaku.segments = newSegments;

      return zaku;
    },

     unarc : function () {
      this.iterate(function (s, index, x, y) {
        var new_segments, nextX, nextY, result = [], name = s[0];

        // Skip anything except arcs
        if (name !== 'A' && name !== 'a') { return null; }

        if (name === 'a') {
          // convert relative arc coordinates to absolute
          nextX = x + s[6];
          nextY = y + s[7];
        } else {
          nextX = s[6];
          nextY = s[7];
        }

        new_segments = a2c(x, y, nextX, nextY, s[4], s[5], s[1], s[2], s[3]);

        // Degenerated arcs can be ignored by renderer, but should not be dropped
        // to avoid collisions with `S A S` and so on. Replace with empty line.
        if (new_segments.length === 0) {
          return [ [ s[0] === 'a' ? 'l' : 'L', s[6], s[7] ] ];
        }

        new_segments.forEach(function (s) {
          result.push([ 'C', s[2], s[3], s[4], s[5], s[6], s[7] ]);
        });

        return result;
      });

      return zaku;
    },

    rel : function () {

      this.iterate(function (s, index, x, y) {
        var name = s[0],
            nameLC = name.toLowerCase(),
            i;

        // Skip relative commands
        if (name === nameLC) { return; }

        // Don't touch the first M to avoid potential confusions.
        if (index === 0 && name === 'M') { return; }

        s[0] = nameLC;

        switch (name) {
          case 'V':
            // V has shifted coords parity
            s[1] -= y;
            return;

          case 'A':
            // ARC is: ['A', rx, ry, x-axis-rotation, large-arc-flag, sweep-flag, x, y]
            // touch x, y only
            s[6] -= x;
            s[7] -= y;
            return;

          default:
            for (i = 1; i < s.length; i++) {
              s[i] -= i % 2 ? x : y; // odd values are X, even - Y
            }
        }
      }, true);

      return zaku;
    }


    } //End of proto

    return this.zakarc.prototype;
 
  },

  _stringify:function(_segments)
    {

        var stringpath = '';
        for (var i = 0; i < _segments.length; i++)
         {
            var _segment = _segments[i];

            for (var j = 0; j < _segment.length; j++)
             {
                stringpath+= (' '+_segment[j]);
             }
         }

       return stringpath;

    },

	_split:function(path)
	 {
      
      /*Convert arc to B-curves*/
        var zakarc = this.zakarc(path)
       
        zakarc.unarc()
        zakarc.rel()
        var segments = zakarc.round(2).segments;
      /*Stringify*/
        var convertedPath = this._stringify(segments);

        var pathSegments = convertedPath.match(pathSegmentPattern);
        return pathSegments;

	},
  split:function(path)
  {
      var newpathString = this.normalizePath(path);
      var newpathArray = this._split(newpathString);
      
      /* Normalize path */
      newpathArray = zaku.normalize(newpathArray);

      return newpathArray;
  },


   ///////////////ANIMATIONS UTILS ////////////////////////


  addkeyframe : function(start, end)
  {
    
    this.utils.keyframe[0] = start;
    this.utils.keyframe[1] = end;

  },

  setInitialValue : function(value)
  {
    
    this.utils.initialValue = value;

  },

  setnewValue : function(value)
  {
    
    this.utils.newValue = value;
    
  },

  setupAnimation(idElem, attr)
  {
      var animation = document.createElementNS('http://www.w3.org/2000/svg' , 'animate');

      animation.setAttribute('href','#'+idElem)
      animation.setAttribute('id','animid')
      animation.setAttribute('attributeName',attr)
      animation.setAttribute('begin','0s')
      animation.setAttribute('from', this.utils.initialValue)
      animation.setAttribute('to', this.utils.newValue)
      animation.setAttribute('dur', (this.utils.keyframe[1] - this.utils.keyframe[0]) + 's')
      animation.setAttribute('fill', 'freeze')
      animation.setAttribute('repeatCount','0')
      
      var svg = document.getElementsByTagName('svg')[0];
      
      svg.appendChild(animation);

  },

   ////////////////////////////////////////////////////////

	place:function(x,y,indexInPathArray, pointIndex, type,fill='tomato')
	{  
    /*make the parallel with c|b and fill style*/
    var natureTemp = fill;
    if(fill === 'b') fill = 'tomato'; 
    if(fill === 'c') fill = 'lightgreen'; 

		var use = document.createElementNS("http://www.w3.org/2000/svg", 'circle');
		use.setAttribute('class','zaku-move-point');
		use.setAttribute('cx', x);
		use.setAttribute('cy', y);
		use.setAttribute('r', '4');
		use.setAttribute('fill', fill);
		use.setAttribute('stroke', '#000');
		use.setAttribute('stroke-width', '2');
    use.setAttribute('type', type);
    use.setAttribute('nature', natureTemp);
    use.setAttribute('pathArrayIndex', indexInPathArray);
		use.setAttribute('pointIndex', pointIndex);

		document.getElementsByTagName('svg')[0].appendChild(use);
		var self = this;
		use.onclick = function(){

			// self.dragstart();
		}
	},
	update: function(elem, pathArray)
	{

       var pathString = pathArray.join(' ');
       elem.setAttribute('d',pathString)

	},
  removepoints:function()
  {
     /**********Clear all movePoints************/
         
       var svg = document.getElementsByTagName('svg')[0];
       var movePoints = document.getElementsByClassName('zaku-move-point');
       var len = movePoints.length
       for (var i = 0; i < len; i++)
        {
          svg.removeChild(movePoints[i])
          //Fixing index issue
          i--;
          len--;
        }

     /*****************************************/

  },
  /*Specifies if H or V*/
  isHCommand:function(firstLetter)
  {
      var _isH = false;

      if(firstLetter === 'h' || firstLetter === 'L')
      {
          _isH = true;
      }
     return _isH;
  },
  isMonoVal:function(firstLetter)
  {
       var _is = false;
      
       switch(firstLetter)
       {

          case 'h' : _is = true; break; case 'H' : _is = true; break; 
          case 'v' : _is = true; break; case 'V' : _is = true; break; 
          default:break;
       }

       return _is;
  },

  isBiVal:function(firstLetter)
  {
     var _is = false;
    
     switch(firstLetter)
     {

        case 'm' : _is = true; break; case 'M' : _is = true; break; 
        case 'l' : _is = true; break; case 'L' : _is = true; break; 
        case 't' : _is = true; break; case 'T' : _is = true; break; 
        default:break;
     }

     return _is;
  },

  isQuadriVal:function(firstLetter)
  {
     var _is = false;
    
     switch(firstLetter)
     {

        case 's' : _is = true; break; case 'S' : _is = true; break; 
        case 'q' : _is = true; break; case 'Q' : _is = true; break; 
        default:break;
     }

     return _is;
  },

  isSixVal:function(firstLetter)
  {
     var _is = false;
    
     switch(firstLetter)
     {

        case 'c' : _is = true; break; case 'C' : _is = true; break; 
        default:break;
     }

     return _is;
  },

  isSevenVal:function(firstLetter)
  {
     var _is = false;
    
     switch(firstLetter)
     {

        case 'a' : _is = true; break; case 'A' : _is = true; break; 
        default:break;
     }

     return _is;
  },

  biValNature:function(firstLetter)
  {
     var _nature = 'b';
    
     switch(firstLetter)
     {

        case 't' : _nature = 'c'; break; case 'T' : _nature = 'c'; break; 
        default:break;
     }

     return _nature;
  },
	setpoints:function(elem, pathArray)
	{
           
        for (var i = 0; i < pathArray.length; i++)
         {
        	
        	var firstLetter = pathArray[i][0];
          // ----------------------------------------//
          if(this.isMonoVal(firstLetter))
          {
              this.placeMonoVal(pathArray , i , firstLetter);
          }

          if(this.isBiVal(firstLetter))
          {
              this.placeBiVal(pathArray , i , firstLetter);
          }

          if(this.isQuadriVal(firstLetter))
          {
              this.placeQuadriVal(pathArray , i , firstLetter);
          }

          if(this.isSixVal(firstLetter))
          {
              this.placeSixVal(pathArray , i , firstLetter);
          }

          // if(this.isSixVal(firstLetter))
          // {
          //     this.placeSixVal(pathArray , i , firstLetter);
          // }


      }

	},

  getprevcoords:function(pathArray, startingIndex)
  {

     /*Attempt to get previous Uppercase command coords*/
     var coordsArray = [];

     for (var i = startingIndex - 1; i >= 0; i--)
      {

         coordsArray.push(pathArray[i])
           /*Stop when meet Majuscule command*/
        if(pathArray[i][0] === pathArray[i][0].toUpperCase())
        {
           // console.log(coordsArray)
           return coordsArray;
        }

        
      }
      
     
  },
  updatepoints:function(elem, pathArray)
  {
           
          this.removepoints();
          
          this.setpoints(elem, pathArray);
          
  },
	 
  setcpointNature:function(_collector,index)
   {
      /* According to modulo we can know if the point is a controlPoint or basicPoint*/
      /* We use '3' (|-->_<--|) because C command contains 3 points*/
       var modulo = index % 3;

       switch(modulo)
        {

           case 0 : _collector['nature'] = 'c'; break;
           case 1 : _collector['nature'] = 'c'; break;
           case 2 : _collector['nature'] = 'b'; break;
           default:break;

        }

       return _collector;

   },
  setQuadriNature:function(_collector,index)
   {
      /* According to modulo we can know if the point is a controlPoint or basicPoint*/
      /* We use '3' (|-->_<--|) because C command contains 3 points*/
       var modulo = index % 2;

       switch(modulo)
        {

           case 0 : _collector['nature'] = 'c'; break;
           case 1 : _collector['nature'] = 'b'; break;
           default:break;

        }

       return _collector;

   },

  /*Rewrite path with for each command the basic number of points */
  /*Transform this "l 25,30 45,80" to this "l 25,30 l 45,80" */
  /*This will simplify relative points calculation : Else (45,80) will use (l 25,30)Not easy to manage */
  normalizePath:function(path)
  {
     var pathArray = this._split(path);
     
     var collectors = [];
     /*Get collection of collectors*/
     for (var i = 0; i < pathArray.length; i++)
      {
         var _command =  pathArray[i][0];
       
         /*Collector for Mono-values commands*/
         if(this.isMonoVal(_command))
         {

            var _collector = this.getMonoValPoint(pathArray,i,_command);

            collectors.push(_collector);

         }

         /*Collector for Bi-values commands*/
         if(this.isBiVal(_command))
         {
            var _collector = this.getBiValPoint(pathArray,i,_command);

            collectors.push(_collector);

         }

         /*Collector for Quadri-values commands*/
         if(this.isQuadriVal(_command))
         {
            var _collector = this.getQuadriValPoint(pathArray,i,_command);

            collectors.push(_collector);

         }


         /*Collector for Six-values commands*/
         if(this.isSixVal(_command))
         {
            var _collector = this.getSixValPoint(pathArray,i,_command);
            // console.log(_collector);
            collectors.push(_collector);

         }



      }

    /*Normalizing all collectors*/
    for (var c = 0; c < collectors.length; c++)
     {
         var collector = collectors[c];
         var command = pathArray[c][0];

         var collectorString = '';
         var space = ' ';
         /*Parsing for Mono*/
         if(this.isMonoVal(command))
         {
            for (var i = 0; i < collector.length; i++) 
             { 

               if (i === collector.length-1) space = '';
               collectorString = collectorString + command + collector[i] + space;
             
             }
         }
        
         /*Parsing for Bi-*/
        if(this.isBiVal(command))
         {
             for (var i = 0; i < collector.length; i++) 
              { 

               if (i === collector.length-1) space = '';

               collectorString = collectorString + command + collector[i][0] + ',' +collector[i][1] + space;
             
              }
          }

         /*Parsing for Bi-*/
        if(this.isQuadriVal(command))
         {   
             var innerspace = ' '; 
             var spaceKiller = 0; 
             for (var i = 0; i < collector.length; i++) 
              { 
               if (i === collector.length-1) innerspace = '';

               collectorString = collectorString + collector[i][0] + ',' +collector[i][1] + innerspace;
             
              }

              collectorString = command+collectorString;

          }

         /*Parsing for Six-*/
        if(this.isSixVal(command))
         {
          for (var m = 0; m < collector.length; m++)
           {
             var innerspace = ' ';
             var _collector = collector[m];
             var innerCollectorString = '';
            
             for (var i = 0; i < _collector.length; i++) 
              { 

                if (i === _collector.length-1) innerspace = '';

                    innerCollectorString = innerCollectorString + _collector[i][0] + ',' +_collector[i][1] + innerspace;
              }
                
              innerCollectorString = command+innerCollectorString;
              
              collectorString+=innerCollectorString;

           }

         }


         /*Update pathArray at index*/
         pathArray[c] = collectorString;


     }
 
     var pathArrayString = pathArray.join('');
     
      
     return pathArrayString;

  },

  nospacepoint:function(coords)
  {
      var len = coords.length;
      for (var i = 0; i < len; i++) 
      {
         if(coords[i] == "")
         {
           coords.splice(i,1);
           i--;
           len--;
         }
      }
      return coords;
  },


  //////////////////////GETTERS //////////////////////////
  /*Mono-value command [V and H]*/

  getMonoValPoint:function(pathArray, indexInPathArray ,command)
  {

          var _coords = [];

          var firstLetter = pathArray[indexInPathArray][0];/*getting command for comparison*/
          var collector = [];

          if(firstLetter === command || firstLetter === command.toUpperCase())
          {

              var coords = pathArray[indexInPathArray].slice(1);

              coords = coords.split(/[, ]/);

              /*Output : [0:"25",
                          1:"75",
                          2:"20",
                          3:"45"]*/
               /*Removing non-desired points (spaces) [..., 4:""]*/
              coords = this.nospacepoint(coords);

              /**************************************/
              /*Get coords couples (x or y): Collecting points(x or y) */
              collector =  coords;
              // console.log(collector)
               

          }

      return collector;
    
  },

  /*Bi-values command  [M , L and T]*/

  getBiValPoint:function(pathArray, indexInPathArray,command)
  {

          var _coords = [];
          

          var firstLetter = pathArray[indexInPathArray][0];/*getting command for comparison*/
          var collector = [];

          if(firstLetter === command || firstLetter === command.toUpperCase())
          {

              var coords = pathArray[indexInPathArray].slice(1);

              coords = coords.split(/[, ]/);
              /*Output : [0:"25",
                          1:"75",
                          2:"20",
                          3:"45"]*/
               /*Removing non-desired points (spaces) [..., 4:""]*/
              coords = this.nospacepoint(coords);
              /**************************************/
              /*Get coords couples (x,y): Collecting points(x,y) */
               var index = 0;
           
               for (var i = 0; i < coords.length; i=i+2)
               {  
                    collector[index] = [];

                     for (var j = i; j < i+2; j++)
                      {

                         collector[index].push(parseFloat(coords[j] , 10));

                      }

                    index++;

               }

            }

      return collector;

  },

  /*Quadri-values command  [S,Q]*/

  getQuadriValPoint:function(pathArray, indexInPathArray,command)
  {

          var _coords = [];
          

          var firstLetter = pathArray[indexInPathArray][0];/*getting command for comparison*/
          var collector = [];

          if(firstLetter === command || firstLetter === command.toUpperCase())
          {

              var coords = pathArray[indexInPathArray].slice(1);

              coords = coords.split(/[, ]/);
              /*Output : [0:"25",
                          1:"75",
                          2:"20",
                          3:"45"]*/
               /*Removing non-desired points (spaces) [..., 4:""]*/
              coords = this.nospacepoint(coords);
              /**************************************/
              /*Get coords couples (x,y): Collecting points(x,y) */
               var index = 0;
           
               for (var i = 0; i < coords.length; i=i+2)
               {  
                    collector[index] = [];
                    collector[index] = this.setQuadriNature(collector[index],index);

                     for (var j = i; j < i+2; j++)
                      {

                         collector[index].push(parseFloat(coords[j] , 10));

                      }

                    index++;

               }

            }

      return collector;

  },

 /*Six-values command  [C]*/

  getSixValPoint:function(pathArray, indexInPathArray,command)
  {

          var _coords = [];
          
          var firstLetter = pathArray[indexInPathArray][0];/*getting command for comparison*/
          var collector = [];
          var collectorArray = [];

          if(firstLetter === command || firstLetter === command.toUpperCase())
          {

              var coords = pathArray[indexInPathArray].slice(1);

              coords = coords.split(/[, ]/);
              /*Output : [0:"25",
                          1:"75",
                          2:"20",
                          3:"45"]*/
               /*Removing non-desired points (spaces) [..., 4:""]*/
              coords = this.nospacepoint(coords);
              // ************************************
              /*Get coords couples (x,y): Collecting points(x,y) */
               var index = 0;
               
               for (var j = 0; j < coords.length; j=j+2)
               {  
                    collector[index] = [];
                    collector[index] = this.setcpointNature(collector[index],index);

                     for (var k = j; k < j+2; k++)
                      {

                         collector[index].push(parseFloat(coords[k] , 10));

                      }

                    index++;

               }
                
              for (var i = 0; i < collector.length/6; i++) 
              {
                 for (j=0 ; j < collector.length; j+=3)
                  {
                      
                      collectorTemp = collector.slice(j,j+3);
                      collectorArray.push(collectorTemp);
                      // do whatever
                  }
              }



        }


      return collectorArray;

  },








  ///////////////// GETTERS END /////////////////////////


  //////////////////////PLACERS //////////////////////////

   placeMonoVal:function(pathArray , indexInPathArray , firstLetter)
   {
      
      var pointsString = pathArray[indexInPathArray].slice(1);
      var coords = [], points = [];
      /*Cpoints : command points*/
      Cpoints = this.getMonoValPoint(pathArray, indexInPathArray,firstLetter.toLowerCase());
      
      for (var j = 0; j < Cpoints.length; j++)
       {
          var pointNature = 'b'; /* By default 'b' */
          var cpoint = Cpoints[j];
        //Setting up the point to place 
         if(this.isHCommand(firstLetter))  
           {
              points[0] = parseFloat(cpoint,10);
              points[1] = 0;
           }
           else
           {
              points[0] = 0;
              points[1] = parseFloat(cpoint,10); 
           }

          /*Get relative coords for positioning in the case of relative command*/
          if (firstLetter == firstLetter.toLowerCase())
            {
             
                 /*Getting previous commands coords*/
                 var _prevcoordsArray = this.getprevcoords(pathArray , indexInPathArray);

                 for (var k = _prevcoordsArray.length - 1; k >= 0; k--)
                  {

                     var _prevcoords = [];
                     var prevcoords = _prevcoordsArray[k].slice(1);

                     var command = _prevcoordsArray[k][0];
                     
                      /*For each command get only the X and Y coord if exist*/
 
                          /*Get MonoVal collector*/

                          if(this.isMonoVal(command))  /*H, V*/
                            {
                                var _prevCollector = this.getMonoValPoint(_prevcoordsArray, k , command);
                                 
                                /*For each point in collector: We normally have just One*/
                                for (var l = 0; l < _prevCollector.length; l++)
                                 { 
                                   /*Update coords : Loop*/
                                   /*To know if we update just X or Y value for H or V*/
                                     if(this.isHCommand(firstLetter)) /*firstLetter = first command*/
                                     {  
                                      
                                       if(this.isHCommand(command)) /**if prevpoint same than dragged point*/
                                        { 
                                          points[0]+= parseFloat(_prevCollector[l], 10);
                                        }
                                        if(!(this.isHCommand(command))) /**if prevpoint same than dragged point*/
                                        { 
                                          points[1]+= parseFloat(_prevCollector[l], 10);
                                        }

                                     }
                                     else
                                     {

                                       if(this.isHCommand(command)) /**if prevpoint same than dragged point*/
                                         { 
                                            points[0]+= parseFloat(_prevCollector[l], 10);
                                         }
                                       if(!(this.isHCommand(command))) /**if prevpoint same than dragged point*/
                                         { 
                                            points[1]+= parseFloat(_prevCollector[l], 10);
                                         }
                                     }
                                 
                                 }
                                
                            }


                          /*Get BiVal collector*/
                          if(this.isBiVal(command))  /*M,L,T*/
                            {
                                var _prevCollector = this.getBiValPoint(_prevcoordsArray, k , command);

                                /*For each point in collector: We normally have just One*/
                                for (var l = 0; l < _prevCollector.length; l++)
                                 { 
                                   /*Update coords : Loop*/
                                   /*To know if we update just X or Y value for H or V*/
                                     if(this.isHCommand(firstLetter)) /*firstLetter = first command*/
                                     {
                                        points[0]+= parseFloat(_prevCollector[l][0], 10);
                                        points[1]+= parseFloat(_prevCollector[l][1], 10);
                                     }
                                     else
                                     {
                                        points[0]+= parseFloat(_prevCollector[l][0], 10);
                                        points[1]+= parseFloat(_prevCollector[l][1], 10);
                                     }
                                 
                                 }
                                
                            }


                          if(this.isSixVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getSixValPoint(_prevcoordsArray, k , command);
                                // console.log(_prevCollector)
                                for (var l = 0; l < _prevCollector.length; l++)
                                 {
                                   /*Next c command points only inherit from the last point of the previous c command :Ref SVG Path c command*/
                                   /*This explains the 2 in _prevCollector[l][2]...*/
                                      points[0]+= parseFloat(_prevCollector[l][2][0], 10);
                                      points[1]+= parseFloat(_prevCollector[l][2][1], 10);

                                  }
                            }

                       
                    }

              }


          this.place(points[0], points[1], indexInPathArray, j, firstLetter, pointNature);

       }

   },
 
   placeBiVal:function(pathArray , indexInPathArray , firstLetter)
   {

      
      var pointsString = pathArray[indexInPathArray].slice(1);
      var coords = [], points = [];
      /*Cpoints : command points*/
      Cpoints = this.getBiValPoint(pathArray, indexInPathArray,firstLetter.toLowerCase());

      for (var j = 0; j < Cpoints.length; j++)
       {
          var pointNature = this.biValNature(firstLetter);
          var cpoint = Cpoints[j];
          
          /*Get relative coords for positioning*/
          if (firstLetter == firstLetter.toLowerCase())
            {
             
                 /*Getting previous commands coords*/
                 var _prevcoordsArray = this.getprevcoords(pathArray , indexInPathArray);
                 // console.log(_prevcoordsArray);
                 for (var k = _prevcoordsArray.length - 1; k >= 0; k--)
                  {

                     var _prevcoords = [];
                     var prevcoords = _prevcoordsArray[k].slice(1);

                      var command = _prevcoordsArray[k][0];
                     
                      /*For each command get only the X and Y coord if exist*/
                           
                           if(this.isMonoVal(command))  /*H, V*/
                            {
                                var _prevCollector = this.getMonoValPoint(_prevcoordsArray, k , command);
                                 
                                /*For each point in collector: We normally have just One*/
                                for (var l = 0; l < _prevCollector.length; l++)
                                 { 
                                   /*Update coords : Loop*/
                                   /*To know if we update just X or Y value for H or V*/
                                     if(this.isHCommand(command)) /*firstLetter = first command*/
                                     {  
                                        cpoint[0]+= parseFloat(_prevCollector[l], 10);
                                     }
                                     else
                                     {
                                        cpoint[1]+= parseFloat(_prevCollector[l], 10);
                                     }
                                 
                                 }
                                
                            }


                          if(this.isBiVal(command))  /*M,L,T*/
                            {
                                var _prevCollector = this.getBiValPoint(_prevcoordsArray, k , command);
                                for (var l = 0; l < _prevCollector.length; l++)
                                 {
                                   /*Update coords : Loop*/
                                    cpoint[0]+= parseFloat(_prevCollector[l][0], 10);
                                    cpoint[1]+= parseFloat(_prevCollector[l][1], 10);
                                 }
                                
                            }
                       
                           if(this.isSixVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getSixValPoint(_prevcoordsArray, k , command);
                                // console.log(_prevCollector)
                                for (var l = 0; l < _prevCollector.length; l++)
                                 {
                                   /*Next c command points only inherit from the last point of the previous c command :Ref SVG Path c command*/
                                   /*This explains the 2 in _prevCollector[l][2]...*/
                                      cpoint[0]+= parseFloat(_prevCollector[l][2][0], 10);
                                      cpoint[1]+= parseFloat(_prevCollector[l][2][1], 10);

                                  }
                            }
                      /*Update coords : Loop*/
                      // cpoint[0]+=_prevcoords[0];
                      // cpoint[1]+=_prevcoords[1];


                 }

              }


          this.place(cpoint[0], cpoint[1], indexInPathArray, j, firstLetter, pointNature);

       }

   },


  placeQuadriVal:function(pathArray , indexInPathArray , firstLetter)
   {

      
      var pointsString = pathArray[indexInPathArray].slice(1);
      var coords = [], points = [];
      /*Cpoints : command points*/
      Cpoints = this.getQuadriValPoint(pathArray, indexInPathArray,firstLetter.toLowerCase());

      for (var j = 0; j < Cpoints.length; j++)
       {
          var cpoint = Cpoints[j];
          
          /*Get relative coords for positioning*/
          if (firstLetter == firstLetter.toLowerCase())
            {
             
                 /*Getting previous commands coords*/
                 var _prevcoordsArray = this.getprevcoords(pathArray , indexInPathArray);
                 // console.log(_prevcoordsArray);
                 for (var k = _prevcoordsArray.length - 1; k >= 0; k--)
                  {

                     var _prevcoords = [];
                     var prevcoords = _prevcoordsArray[k].slice(1);

                      var command = _prevcoordsArray[k][0];
                     
                      /*For each command get only the X and Y coord if exist*/
                           
                           if(this.isMonoVal(command))  /*H, V*/
                            {
                                var _prevCollector = this.getMonoValPoint(_prevcoordsArray, k , command);
                                 
                                /*For each point in collector: We normally have just One*/
                                for (var l = 0; l < _prevCollector.length; l++)
                                 { 
                                   /*Update coords : Loop*/
                                   /*To know if we update just X or Y value for H or V*/
                                     if(this.isHCommand(command)) /*firstLetter = first command*/
                                     {  
                                        cpoint[0]+= parseFloat(_prevCollector[l], 10);
                                     }
                                     else
                                     {
                                        cpoint[1]+= parseFloat(_prevCollector[l], 10);
                                     }
                                 
                                 }
                                
                            }


                          if(this.isBiVal(command))  /*M,L,T*/
                            {
                                var _prevCollector = this.getBiValPoint(_prevcoordsArray, k , command);
                                for (var l = 0; l < _prevCollector.length; l++)
                                 {
                                   /*Update coords : Loop*/
                                    cpoint[0]+= parseFloat(_prevCollector[l][0], 10);
                                    cpoint[1]+= parseFloat(_prevCollector[l][1], 10);
                                 }
                                
                            }


                           if(this.isQuadriVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getQuadriValPoint(_prevcoordsArray, k , command);
                                 
                                      /*Next command points only inherit from the last point of the previous s command :Ref SVG Path c command*/
                                      /*This explains the 2 in _prevCollector[l][2]...*/
                                      cpoint[0]+= parseFloat(_prevCollector[1][0], 10);
                                      cpoint[1]+= parseFloat(_prevCollector[1][1], 10);

                            }

                       
                           if(this.isSixVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getSixValPoint(_prevcoordsArray, k , command);
                                // console.log(_prevCollector)
                                for (var l = 0; l < _prevCollector.length; l++)
                                 {
                                   /*Next c command points only inherit from the last point of the previous c command :Ref SVG Path c command*/
                                   /*This explains the 2 in _prevCollector[l][2]...*/
                                      cpoint[0]+= parseFloat(_prevCollector[l][2][0], 10);
                                      cpoint[1]+= parseFloat(_prevCollector[l][2][1], 10);

                                  }
                            }
                   

                  }

              }

            this.place(cpoint[0], cpoint[1], indexInPathArray, j, firstLetter, cpoint.nature);

       }

   },



  placeSixVal:function(pathArray , indexInPathArray , firstLetter)
   {

      var pointsString = pathArray[indexInPathArray].slice(1);
      var coords = [], points = [];
      /*Cpoints : command points*/
      Cpoints = this.getSixValPoint(pathArray, indexInPathArray,firstLetter.toLowerCase());

      for (var j = 0; j < Cpoints.length; j++)
       {
          // var pointNature = this.biValNature(firstLetter);
          var pointNature = 'c';
          var cpoint = Cpoints[j];

          // console.log(cpoint);

         for (var i = 0; i < cpoint.length; i++) 
         {
          
           var innerCpoint = cpoint[i];
           
           // console.log(innerCpoint)
          
          /*Get relative coords for positioning*/
           if (firstLetter == firstLetter.toLowerCase())
             {
               

                 /*Getting previous commands coords*/
                 var _prevcoordsArray = this.getprevcoords(pathArray , indexInPathArray);
                 // console.log(_prevcoordsArray);
                 for (var k = _prevcoordsArray.length - 1; k >= 0; k--)
                  {

                     var _prevcoords = [];
                     var prevcoords = _prevcoordsArray[k].slice(1);
                      

                      var command = _prevcoordsArray[k][0];
                      /*For each command get only the X and Y coord if exist*/
                           
                           if(this.isMonoVal(command))  /*H, V*/
                            {
                                var _prevCollector = this.getMonoValPoint(_prevcoordsArray, k , command);
                                 
                                /*For each point in collector: We normally have just One*/
                                for (var l = 0; l < _prevCollector.length; l++)
                                 { 
                                   /*Update coords : Loop*/
                                   /*To know if we update just X or Y value for H or V*/
                                     if(this.isHCommand(command)) //firstLetter = first command
                                     {  
                                        innerCpoint[0]+= parseFloat(_prevCollector[l], 10);
                                     }
                                     else
                                     {
                                        innerCpoint[1]+= parseFloat(_prevCollector[l], 10);
                                     }
                                 
                                 }
                                
                            }


                          if(this.isBiVal(command))  /*M,L,T*/
                            {
                                var _prevCollector = this.getBiValPoint(_prevcoordsArray, k , command);
                                for (var l = 0; l < _prevCollector.length; l++)
                                 {
                                   /*Update coords : Loop*/
                                    innerCpoint[0]+= parseFloat(_prevCollector[l][0], 10);
                                    innerCpoint[1]+= parseFloat(_prevCollector[l][1], 10);
                                 }
                                
                            }


                           if(this.isQuadriVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getQuadriValPoint(_prevcoordsArray, k , command);
                               
                                    /*Next command points only inherit from the last point of the previous s command :Ref SVG Path c command*/
                                    /*This explains the 2 in _prevCollector[l][2]...*/
                                    innerCpoint[0]+= parseFloat(_prevCollector[1][0], 10);
                                    innerCpoint[1]+= parseFloat(_prevCollector[1][1], 10);

                            }

                          if(this.isSixVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getSixValPoint(_prevcoordsArray, k , command);
                                // console.log(_prevCollector)
                                for (var l = 0; l < _prevCollector.length; l++)
                                 {
                                   /*Next c command points only inherit from the last point of the previous c command :Ref SVG Path c command*/
                                   /*This explains the 2 in _prevCollector[l][2]...*/
                                      innerCpoint[0]+= parseFloat(_prevCollector[l][2][0], 10);
                                      innerCpoint[1]+= parseFloat(_prevCollector[l][2][1], 10);

                                  }
                            }

                     

                 }

              } /*If relative command end*/

          // //Check index in pathArray
          this.place(innerCpoint[0], innerCpoint[1], indexInPathArray, i , firstLetter, innerCpoint.nature);

         }//For each point in C command array

       }

   },


  //////////////////////PLACERS END //////////////////////////


  ///////////////// SETTERS  ////////////////////////////

setMonoValPoint:function(pathArray, pathArrayIndex, pointIndex, type ,mpoints, newX, newY, movePointX, movePointY)
   {
      // console.log(mpoints)
      var firstLetter = pathArray[pathArrayIndex][0];
      /*Updating Point*/
      //Setting up the point to update 
         if(this.isHCommand(firstLetter))  
           {
              mpoints[pointIndex] = newX+movePointX;
           }
           else
           {

              mpoints[pointIndex] = newY+movePointY;
           }
          
     
      /*Updating mpoints Array*/
      /*Updating pathArray*/
      var pointsString = '';
      var space = ' ';
      /*Normally just one point with one coord (X or Y)*/
      for (var j = 0; j < mpoints.length; j++)
       {
        
        /*Get relative coords for positioning for the Dragged Point */ 
         if((firstLetter == firstLetter.toLowerCase()) && (j == pointIndex))
            {
             
                 /*Getting previous commands coords*/
                 var _prevcoordsArray = this.getprevcoords(pathArray , pathArrayIndex);

                 for (var k = _prevcoordsArray.length - 1; k >= 0; k--)
                  {

                     var _prevcoords = [];
                     var prevcoords = _prevcoordsArray[k].slice(1);

                      var command = _prevcoordsArray[k][0];
                     
                      /*For each command get only the X and Y coord of its points if exist*/

                      /*Get MonoVal collector*/
                      if(this.isMonoVal(command))  /*H, V*/
                        {
                            var _prevCollector = this.getMonoValPoint(_prevcoordsArray, k , command);
                             
                            /*For each point in collector: We normally have just One*/
                            for (var l = 0; l < _prevCollector.length; l++)
                             { 
                               /*Update coords : Loop*/
                               /*To know if we update just X or Y value for H or V*/
                                 if(this.isHCommand(firstLetter)) /*firstLetter = first command*/
                                 {  
                                    if(this.isHCommand(command)) /*if prevpoint same than dragged point*/
                                    { /*V have not effect on H coords calculation : vice versa */
                                       mpoints[j]-= parseFloat(_prevCollector[l], 10);
                                    }
                                 }
                                 else
                                 {
                                   if(!(this.isHCommand(command))) /*if prevpoint same than dragged point*/
                                    { 
                                       mpoints[j]-= parseFloat(_prevCollector[l], 10);
                                    }
                                 }
                             
                             }
                            
                        }


                        if(this.isBiVal(command))  /*M,L,T*/
                        {
                            var _prevCollector = this.getBiValPoint(_prevcoordsArray, k , command);
                            // console.log(_prevCollector)
                            for (var l = 0; l < _prevCollector.length; l++)
                             { 
                                 /*Update coords : Loop*/ 
                                 /*To know if we update just X or Y value for H or V*/
                                   if(this.isHCommand(firstLetter)) /*firstLetter = first command*/
                                   {  
                                      mpoints[j]-= parseFloat(_prevCollector[l][0], 10);
                                   }
                                   else
                                   {
                                      mpoints[j]-= parseFloat(_prevCollector[l][1], 10);
                                   }
                                 
                             }
                            
                        }


                       if(this.isQuadriVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getQuadriValPoint(_prevcoordsArray, k , command);
                                 
                                  /*Next command points only inherit from the last point of the previous s command :Ref SVG Path c command*/
                                  /*This explains the 2 in _prevCollector[l][2]...*/
                                  
                                  if(this.isHCommand(firstLetter)) /*firstLetter = first command*/
                                   {  
                                      mpoints[j]-= parseFloat(_prevCollector[1][0], 10);
                                   }
                                   else
                                   {
                                      mpoints[j]-= parseFloat(_prevCollector[1][1], 10);
                                   }
                               

                            }

 
                        if(this.isSixVal(command))  /*C*/
                          {   
                              var _prevCollector = this.getSixValPoint(_prevcoordsArray, k , command);
                              // console.log(_prevCollector)
                              for (var l = 0; l < _prevCollector.length; l++)
                               {
                                 /*Next c command points only inherit from the last point of the previous c command :Ref SVG Path c command*/
                                 /*This explains the 2 in _prevCollector[l][2]...*/
                                 /*To know if we update just X or Y value for H or V*/
                                   if(this.isHCommand(firstLetter)) /*firstLetter = first command*/
                                   {  
                                      mpoints[j]-= parseFloat(_prevCollector[l][2][0], 10);
                                   }
                                   else
                                   {
                                      mpoints[j]-= parseFloat(_prevCollector[l][2][1], 10);
                                   }
                                   
                                }

                          }

                     

                 }

            }


          if (j === mpoints.length-1)
           space = ''; /*Avoid adding space at the end of the string : The reverse caused NAN,NAN as new point*/
           pointsString = pointsString + mpoints[j] + space;
       }

      pathArray[pathArrayIndex] = type+pointsString;

      return pathArray;

  },


/*setBiValPoint

   

*/
   
setBiValPoint:function(pathArray, pathArrayIndex, pointIndex, type ,mpoints, newX, newY, movePointX, movePointY)
   {
      // console.log(mpoints)
      var firstLetter = pathArray[pathArrayIndex][0];
      /*Updating Point*/
      mpoints[pointIndex][0] = newX+movePointX;
      mpoints[pointIndex][1] = newY+movePointY;

      /*Updating mpoints Array*/
      /*Updating pathArray*/
      var pointsString = '';
      var space = ' ';
      
      for (var j = 0; j < mpoints.length; j++)
       {
        
        /*Get relative coords for positioning for the Dragged Point among those of C command*/
         if((firstLetter == firstLetter.toLowerCase()) && (j == pointIndex))
            {
             
                 /*Getting previous commands coords*/
                 var _prevcoordsArray = this.getprevcoords(pathArray , pathArrayIndex);

                 for (var k = _prevcoordsArray.length - 1; k >= 0; k--)
                  {

                     var _prevcoords = [];
                     var prevcoords = _prevcoordsArray[k].slice(1);

                      var command = _prevcoordsArray[k][0];
                     
                      /*For each command get only the X and Y coord of its points if exist*/

                        /*Get MonoVal collector*/
                        if(this.isMonoVal(command))  /*H, V*/
                          {
                              var _prevCollector = this.getMonoValPoint(_prevcoordsArray, k , command);
                               
                              /*For each point in collector: We normally have just One*/
                              for (var l = 0; l < _prevCollector.length; l++)
                               { 
                                 /*Update coords : Loop*/
                                 /*To know if we update just X or Y value for H or V*/
                                   if(this.isHCommand(command)) /*firstLetter = first command*/
                                   {  
                                       mpoints[j][0]-= parseFloat(_prevCollector[l], 10);
                                    
                                   }
                                   else
                                   {
                                     
                                      mpoints[j][1]-= parseFloat(_prevCollector[l], 10);
                                      
                                   }
                               
                               }
                              
                          }

                        if(this.isBiVal(command))  /*M,L,T*/
                        {
                            var _prevCollector = this.getBiValPoint(_prevcoordsArray, k , command);
                            // console.log(_prevCollector)
                            for (var l = 0; l < _prevCollector.length; l++)
                             { 
                               /*Update coords : Loop*/
                                mpoints[j][0]-= parseFloat(_prevCollector[l][0], 10);
                                mpoints[j][1]-= parseFloat(_prevCollector[l][1], 10);
                             }
                            
                        }

                        if(this.isQuadriVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getQuadriValPoint(_prevcoordsArray, k , command);
                           
                                /*Next command points only inherit from the last point of the previous s command :Ref SVG Path c command*/
                                /*This explains the 2 in _prevCollector[l][2]...*/
                                mpoints[j][0]-= parseFloat(_prevCollector[1][0], 10);
                                mpoints[j][1]-= parseFloat(_prevCollector[1][1], 10);

                            }


                        if(this.isSixVal(command))  /*C*/
                          {   
                              var _prevCollector = this.getSixValPoint(_prevcoordsArray, k , command);
                              // console.log(_prevCollector)
                              for (var l = 0; l < _prevCollector.length; l++)
                               {
                                 /*Next c command points only inherit from the last point of the previous c command :Ref SVG Path c command*/
                                 /*This explains the 2 in _prevCollector[l][2]...*/
                                    mpoints[j][0]-= parseFloat(_prevCollector[l][2][0], 10);
                                    mpoints[j][1]-= parseFloat(_prevCollector[l][2][1], 10);

                                }

                          }


                       
                      

                 }

            }


          if (j === mpoints.length-1)
           space = ''; /*Avoid adding space at the end of the string : The reverse caused NAN,NAN as new point*/
           pointsString = pointsString + mpoints[j][0] +','+ mpoints[j][1] + space;
       }

      pathArray[pathArrayIndex] = type+pointsString;

      return pathArray;

  },



   
setQuadriValPoint:function(pathArray, pathArrayIndex, pointIndex, type ,mpoints, newX, newY, movePointX, movePointY)
   {
        // console.log(mpoints)
      var firstLetter = pathArray[pathArrayIndex][0];
      /*Updating Point*/
      // console.log(mpoints[pointIndex])
      
      /*mpoints[0] due to array structuration*/
      mpoints[pointIndex][0] = newX+movePointX;
      mpoints[pointIndex][1] = newY+movePointY;
      

      /*Updating mpoints Array*/
      /*Updating pathArray*/
      var pointsString = '';
      var space = ' ';
       
        
        /*Get relative coords for positioning for the Dragged Point among those of C command*/
         if(firstLetter == firstLetter.toLowerCase() )
            {
             
                 /*Getting previous commands coords*/
                 var _prevcoordsArray = this.getprevcoords(pathArray , pathArrayIndex);

                 for (var k = _prevcoordsArray.length - 1; k >= 0; k--)
                  {

                     var _prevcoords = [];
                     var prevcoords = _prevcoordsArray[k].slice(1);

                      var command = _prevcoordsArray[k][0];
                     
                      /*For each command get only the X and Y coord of its points if exist*/

                        /*Get MonoVal collector*/
                        if(this.isMonoVal(command))  /*H, V*/
                          {
                              var _prevCollector = this.getMonoValPoint(_prevcoordsArray, k , command);
                               
                              /*For each point in collector: We normally have just One*/
                              for (var l = 0; l < _prevCollector.length; l++)
                               { 
                                 /*Update coords : Loop*/
                                 /*To know if we update just X or Y value for H or V*/
                                   if(this.isHCommand(command)) /*firstLetter = first command*/
                                   {  
                                       mpoints[pointIndex][0]-= parseFloat(_prevCollector[l], 10);
                                    
                                   }
                                   else
                                   {
                                     
                                      mpoints[pointIndex][1]-= parseFloat(_prevCollector[l], 10);
                                      
                                   }
                               
                               }
                              
                          }

                        if(this.isBiVal(command))  /*M,L,T*/
                        {
                            var _prevCollector = this.getBiValPoint(_prevcoordsArray, k , command);
                            // console.log(_prevCollector)
                            for (var l = 0; l < _prevCollector.length; l++)
                             { 
                               /*Update coords : Loop*/
                                mpoints[pointIndex][0]-= parseFloat(_prevCollector[l][0], 10);
                                mpoints[pointIndex][1]-= parseFloat(_prevCollector[l][1], 10);
                             }
                            
                        }


                        if(this.isQuadriVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getQuadriValPoint(_prevcoordsArray, k , command);
                           
                                /*Next command points only inherit from the last point of the previous s command :Ref SVG Path c command*/
                                /*This explains the 2 in _prevCollector[l][2]...*/
                                mpoints[pointIndex][0]-= parseFloat(_prevCollector[1][0], 10);
                                mpoints[pointIndex][1]-= parseFloat(_prevCollector[1][1], 10);

                            }



                       if(this.isSixVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getSixValPoint(_prevcoordsArray, k , command);
                                // console.log(_prevCollector)
                                for (var l = 0; l < _prevCollector.length; l++)
                                 { 
                                   /*Next c command points only inherit from the last point of the previous c command :Ref SVG Path c command*/
                                   /*This explains the 2 in _prevCollector[l][2]...*/
                                      mpoints[pointIndex][0]-= parseFloat(_prevCollector[l][2][0], 10);
                                      mpoints[pointIndex][1]-= parseFloat(_prevCollector[l][2][1], 10);

                                  }
                            }
                       
                      

                 }

            }
         
 
       //Converting mpoints Array to String 

       for (var m = 0; m < mpoints.length; m++)
           {
             var innerspace = ' ';
             var _collector = mpoints[m];
             
                if (m === _collector.length-1) innerspace = '';

                    pointsString = pointsString + _collector[0] + ',' +_collector[1] + innerspace;
             
           }

          pointsString = type+pointsString;
          // console.log(pointsString)

          pathArray[pathArrayIndex] = pointsString;


      return pathArray;


  },



setSixValPoint:function(pathArray, pathArrayIndex, pointIndex, type ,mpoints, newX, newY, movePointX, movePointY)
   {
      // console.log(mpoints)
      var firstLetter = pathArray[pathArrayIndex][0];
      /*Updating Point*/
      // console.log(mpoints[0][pointIndex])
      
      /*mpoints[0] due to array structuration*/
      mpoints[0][pointIndex][0] = newX+movePointX;
      mpoints[0][pointIndex][1] = newY+movePointY;
      

      /*Updating mpoints Array*/
      /*Updating pathArray*/
      var pointsString = '';
      var space = ' ';
       
        
        /*Get relative coords for positioning for the Dragged Point among those of C command*/
         if(firstLetter == firstLetter.toLowerCase() )
            {
             
                 /*Getting previous commands coords*/
                 var _prevcoordsArray = this.getprevcoords(pathArray , pathArrayIndex);

                 for (var k = _prevcoordsArray.length - 1; k >= 0; k--)
                  {

                     var _prevcoords = [];
                     var prevcoords = _prevcoordsArray[k].slice(1);

                      var command = _prevcoordsArray[k][0];
                     
                      /*For each command get only the X and Y coord of its points if exist*/

                        /*Get MonoVal collector*/
                        if(this.isMonoVal(command))  /*H, V*/
                          {
                              var _prevCollector = this.getMonoValPoint(_prevcoordsArray, k , command);
                               
                              /*For each point in collector: We normally have just One*/
                              for (var l = 0; l < _prevCollector.length; l++)
                               { 
                                 /*Update coords : Loop*/
                                 /*To know if we update just X or Y value for H or V*/
                                   if(this.isHCommand(command)) /*firstLetter = first command*/
                                   {  
                                       mpoints[0][pointIndex][0]-= parseFloat(_prevCollector[l], 10);
                                    
                                   }
                                   else
                                   {
                                     
                                      mpoints[0][pointIndex][1]-= parseFloat(_prevCollector[l], 10);
                                      
                                   }
                               
                               }
                              
                          }

                        if(this.isBiVal(command))  /*M,L,T*/
                        {
                            var _prevCollector = this.getBiValPoint(_prevcoordsArray, k , command);
                            // console.log(_prevCollector)
                            for (var l = 0; l < _prevCollector.length; l++)
                             { 
                               /*Update coords : Loop*/
                                mpoints[0][pointIndex][0]-= parseFloat(_prevCollector[l][0], 10);
                                mpoints[0][pointIndex][1]-= parseFloat(_prevCollector[l][1], 10);
                             }
                            
                        }



                        if(this.isQuadriVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getQuadriValPoint(_prevcoordsArray, k , command);
                           
                                /*Next command points only inherit from the last point of the previous s command :Ref SVG Path c command*/
                                /*This explains the 2 in _prevCollector[l][2]...*/
                                mpoints[j][0]-= parseFloat(_prevCollector[1][0], 10);
                                mpoints[j][1]-= parseFloat(_prevCollector[1][1], 10);

                            }


                       if(this.isSixVal(command))  /*C*/
                            {   
                                var _prevCollector = this.getSixValPoint(_prevcoordsArray, k , command);
                                // console.log(_prevCollector)
                                for (var l = 0; l < _prevCollector.length; l++)
                                 { 
                                   /*Next c command points only inherit from the last point of the previous c command :Ref SVG Path c command*/
                                   /*This explains the 2 in _prevCollector[l][2]...*/
                                      mpoints[0][pointIndex][0]-= parseFloat(_prevCollector[l][2][0], 10);
                                      mpoints[0][pointIndex][1]-= parseFloat(_prevCollector[l][2][1], 10);

                                  }
                            }
                       
                      

                 }

            }

 
       //Converting mpoints Array to String 

       for (var m = 0; m < mpoints.length; m++)
           {
             var innerspace = ' ';
             var _collector = mpoints[m];
             var innerCollectorString = '';
            
             for (var i = 0; i < _collector.length; i++) 
              { 

                if (i === _collector.length-1) innerspace = '';

                    innerCollectorString = innerCollectorString + _collector[i][0] + ',' +_collector[i][1] + innerspace;
              }
                
              innerCollectorString = type+innerCollectorString;
              pointsString = innerCollectorString;

           }

         pathArray[pathArrayIndex] = pointsString;


      return pathArray;

  },

  ///////////////// SETTERS END /////////////////////////

  ///////////////// TESTERS /////////////////////////

  testMonoval:function(movePointX, movePointY)
  {
      var pathArrayIndex = 3;
      var pointIndex = 0;
      var newX = 5;
      var newY = 10;
      var command = 'v';

      var mpoints = this.getMonoValPoint(pathArray, pathArrayIndex , command);

      console.log(mpoints);
      // console.log(pathArray);

      /*Updating point and pathArray*/
      var newpathArray = this.setMonoValPoint(pathArray, pathArrayIndex, pointIndex, command ,mpoints, newX, newY, movePointX, movePointY)
      console.log(newpathArray)
      // var newmpoints = this.setcpoints(pathArray, pathArrayIndex, pointIndex, mpoints,newX, newY, movePointX, movePointY, len);
      // var pointsString = 'C';
      
      

  },

  testBival:function(movePointX, movePointY)
  {
      var pathArrayIndex = 1;
      var pointIndex = 0;
      var newX = 5;
      var newY = 10;
      var command = 'l';

      var mpoints = this.getBiValPoint(pathArray, pathArrayIndex , command);

      console.log(mpoints);
      // console.log(pathArray);

      /*Updating point and pathArray*/
      var newpathArray = this.setBiValPoint(pathArray, pathArrayIndex, pointIndex, command ,mpoints, newX, newY, movePointX, movePointY)
      console.log(newpathArray)
      // var newmpoints = this.setcpoints(pathArray, pathArrayIndex, pointIndex, mpoints,newX, newY, movePointX, movePointY, len);
      // var pointsString = 'C';
      
      

  },

  testQuadrival:function(movePointX, movePointY)
  {
      var pathArrayIndex = 2;
      var pointIndex = 0;
      var newX = 5;
      var newY = 10;
      var command = 'S';

      var mpoints = this.getQuadriValPoint(pathArray, pathArrayIndex , command);

      console.log(mpoints);
      // console.log(pathArray);

      /*Updating point and pathArray*/
      var newpathArray = this.setQuadriValPoint(pathArray, pathArrayIndex, pointIndex, command ,mpoints, newX, newY, movePointX, movePointY)
      console.log(newpathArray)
      // var newmpoints = this.setcpoints(pathArray, pathArrayIndex, pointIndex, mpoints,newX, newY, movePointX, movePointY, len);
      // var pointsString = 'C';
      
      

  },

  testSixval:function(movePointX, movePointY)
  {
      var pathArrayIndex = 6;
      var pointIndex = 0;
      var newX = 5;
      var newY = 10;
      var command = 'c';
      
      var mpoints = this.getSixValPoint(pathArray, pathArrayIndex , command);

      console.log(mpoints);
      // console.log(pathArray);

      /*Updating point and pathArray*/
      var newpathArray = this.setSixValPoint(pathArray, pathArrayIndex, pointIndex, command ,mpoints, newX, newY, movePointX, movePointY)
      console.log(newpathArray)
      // var newmpoints = this.setcpoints(pathArray, pathArrayIndex, pointIndex, mpoints,newX, newY, movePointX, movePointY, len);
      // var pointsString = 'C';
      
      

  },


  normalizePath_test:function()
  {
     
    var newpathArray = this.normalizePath(pathArray);
    console.log(newpathArray)
  },

  startAnim : function(elem,start,end)
  {
      
      /*Keyframe*/
      this.addkeyframe(start , end);
      /*Setting attr init value*/
      var initval = elem.getAttribute('d');
      this.setInitialValue(initval);

  },
  endAnim:function(elem)
  {
       /*Setting attr new value*/
        var newval = elem.getAttribute('d');
        this.setnewValue(newval);
      /*Add animation : */
        var idElem = elem.getAttribute('id');
        this.setupAnimation(idElem, 'd');

  }

}


 
 


