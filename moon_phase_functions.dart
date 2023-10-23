// ignore_for_file: prefer_typing_uninitialized_variables, non_constant_identifier_names, unused_local_variable, duplicate_ignore

import 'dart:math';

import 'package:flutter/material.dart';
import 'package:intl/intl.dart';

import 'common_functions.dart';

class SunMoonTime {
  DateTime date;
  double tz;
  double longitude;
  double latitude;

  SunMoonTime({
    required this.date,
    required this.tz,
    required this.latitude,
    required this.longitude,
  });

  Map moonRiseSet() {
    double mjds = mjd();
    // ignore: unused_local_variable
    var sglong, sglat, date, ym, yz, above, utrise, utset, j;
    var yp, nz, rise, sett, hour, z1, z2, iobj, rads = 0.0174532925;
    List quadout = [];
    var sinho;
    Map outstring = {"rise": false, "set": false, "transit": false};
    sinho = sin(rads * 8 / 60);
    sglat = sin(rads * latitude);
    double cglat = cos(rads * latitude);
    date = mjds - tz / 24;
    rise = false;
    sett = false;
    above = false;
    hour = 1.0;
    ym = sin_alt(1, date, hour - 1.0, longitude, cglat, sglat) - sinho;
    if (ym > 0.0) above = true;
    while (hour < 25 && (sett == false || rise == false)) {
      yz = sin_alt(1, date, hour, longitude, cglat, sglat) - sinho;
      yp = sin_alt(1, date, hour + 1.0, longitude, cglat, sglat) - sinho;
      quadout = quad(ym, yz, yp);
      nz = quadout[0];
      z1 = quadout[1];
      z2 = quadout[2];
      double xe = quadout[3];
      double ye = quadout[4];

      // case when one event is found in the interval
      if (nz == 1) {
        if (ym < 0.0) {
          utrise = hour + z1;
          rise = true;
        } else {
          utset = hour + z1;
          sett = true;
        }
      }

      if (nz == 2) {
        if (ye < 0.0) {
          utrise = hour + z2;
          utset = hour + z1;
        } else {
          utrise = hour + z1;
          utset = hour + z2;
        }
      }

      ym = yp;
      hour += 2.0;
    }

    if (rise == true || sett == true) {
      if (rise == true) {
        String hrsmins = hrsmin(utrise);
        Map time = hrsStringToObject(hrsmins);

        Map a = {
          "hour": time["hour"],
          "minute": time["minute"],
          "meridian": time["meridian"],
          "hrsmin": time["hourMins"],
          "milliSecond": time["milliSecond"],
          "object": time["object"],
        };

        outstring["rise"] = a;
      }
      if (sett == true) {
        String hrsmins = hrsmin(utset);
        Map time = hrsStringToObject(hrsmins);

        Map b = {
          "hour": time["hour"],
          "minute": time["minute"],
          "meridian": time["meridian"],
          "hrsmin": time["hourMins"],
          "milliSecond": time["milliSecond"],
          "object": time["object"],
        };
        outstring["set"] = b;
      }
    }

    ///Moon:
    if (rise != false && sett != false) {
      var risehrsmin = outstring["rise"]["hrsmin"].split(":");
      var sethrsmin = outstring["set"]["hrsmin"].split(":");

      final rise = DateTime(this.date.year, this.date.month, this.date.day, int.parse(risehrsmin[0]), int.parse(risehrsmin[1]));
      final set = DateTime(this.date.year, this.date.month, this.date.day, int.parse(sethrsmin[0]) + 4, int.parse(sethrsmin[1]));
      final transit = rise.add(Duration(milliseconds: (set.millisecondsSinceEpoch - rise.millisecondsSinceEpoch) ~/ 2));

      DateTime time = DateTime(this.date.year, this.date.month, this.date.day, transit.hour, transit.minute, 00);
      Map c = {
        "hour": transit.hour > 12
            ? transit.hour - 12
            : transit.hour == 0
                ? 12
                : transit.hour,
        "minute": transit.minute,
        "meridian": DateFormat("a").format(transit),
        "hrsmin": "${transit.hour}:${transit.minute}",
        "milliSecond": time.millisecondsSinceEpoch,
        "object": time
      };
      outstring["transit"] = c;
    }

    return outstring;
  }

  sunRiseSet() {
    double mjds = mjd();
    double cglat, xe, ye;
    var sglat, date, ym, yz, above, utrise, utset;
    var yp, nz, rise, sett, hour, z1, z2, iobj, rads = 0.0174532925;
    List quadout = [];
    double sinho;
    Map outstring = {"rise": false, "set": false, "transit": false};

    sinho = sin(rads * -0.833);
    sglat = sin(rads * latitude);
    cglat = cos(rads * latitude);
    date = mjds - tz / 24;
    rise = false;
    sett = false;
    above = false;
    hour = 1.0;
    ym = sin_alt(2, date, hour - 1.0, longitude, cglat, sglat) - sinho;
    if (ym > 0.0) above = true;

    while (hour < 25 && (sett == false || rise == false)) {
      yz = sin_alt(2, date, hour, longitude, cglat, sglat) - sinho;
      yp = sin_alt(2, date, hour + 1.0, longitude, cglat, sglat) - sinho;
      quadout = quad(ym, yz, yp);
      nz = quadout[0];
      z1 = quadout[1];
      z2 = quadout[2];
      xe = quadout[3];
      ye = quadout[4];

      if (nz == 1) {
        if (ym < 0.0) {
          utrise = hour + z1;
          rise = true;
        } else {
          utset = hour + z1;
          sett = true;
        }
      }
      if (nz == 2) {
        if (ye < 0.0) {
          utrise = hour + z2;
          utset = hour + z1;
        } else {
          utrise = hour + z1;
          utset = hour + z2;
        }
      }
      ym = yp;
      hour += 2.0;
    }

    if (rise == true || sett == true) {
      if (rise == true) {
        String hrsmins = hrsmin(utrise);
        Map time = hrsStringToObject(hrsmins);

        Map a = {
          "hour": time["hour"],
          "minute": time["minute"],
          "meridian": time["meridian"],
          "hrsmin": time["hourMins"],
          "milliSecond": time["milliSecond"],
          "object": time["object"],
        };

        outstring["rise"] = a;
      }
      if (sett == true) {
        String hrsmins = hrsmin(utset);
        Map time = hrsStringToObject(hrsmins);

        Map b = {
          "hour": time["hour"],
          "minute": time["minute"],
          "meridian": time["meridian"],
          "hrsmin": time["hourMins"],
          "milliSecond": time["milliSecond"],
          "object": time["object"],
        };
        outstring["set"] = b;
      }
    }

    if (rise != false && sett != false) {
      var risehrsmin = outstring["rise"]["hrsmin"].split(":");
      var sethrsmin = outstring["set"]["hrsmin"].split(":");
      var x = int.parse(risehrsmin[0]) * 60 + int.parse(risehrsmin[1]);
      var y = int.parse(sethrsmin[0]) * 60 + int.parse(sethrsmin[1]);

      var z = ((x + y) / 2).floor();
      var d = (z / 60).floor().toString() + (z % 60).toString();
      Map time = hrsStringToObject(d.toString());
      Map c = {
        "hour": time["hour"],
        "minute": time["minute"],
        "meridian": time["meridian"],
        "hrsmin": time["hourMins"],
        "milliSecond": time["milliSecond"],
        "object": time["object"],
      };

      outstring["transit"] = c;
    }

    return outstring;
  }

  hrsStringToObject(String str) {
    Map timeData = {"hour": "00", "minute": "00", "meridian": "AM", "hourMins": "00:00", "milliSecond": "0000000000000", "object": null};

    int position = str.length - 2;

    String timeString = "${str.substring(0, position)}:${str.substring(position)}";
    final dateList = timeString.split(":");
    timeData["hourMins"] = timeString;
    timeData["hour"] = (int.parse(dateList[0]) > 12)
        ? (int.parse(dateList[0]) - 12).toString()
        : (int.parse(dateList[0]) == 0)
            ? "12"
            : dateList[0];
    timeData["minute"] = dateList[1];
    timeData["meridian"] = int.parse(dateList[0]) >= 12 ? "PM" : "AM";

    DateTime time = DateTime(date.year, date.month, date.day, int.parse(dateList[0]), int.parse(dateList[1]), 0);

    timeData["milliSecond"] = time.millisecondsSinceEpoch;
    timeData['object'] = time;

    return timeData;
  }

  dateTimeObjectToHrsObject(DateTime dateTime) {
    Map timeData = {"hour": "00", "minute": "00", "meridian": "AM", "hourMins": "00:00", "milliSecond": "00", "object": null};

    timeData["hourMins"] = "${dateTime.hour}:${dateTime.minute}";
    timeData["hour"] = (dateTime.hour > 12)
        ? (dateTime.hour - 12).toString()
        : (dateTime.hour == 0)
            ? "12"
            : dateTime.hour;
    timeData["minute"] = dateTime.minute;
    timeData["meridian"] = dateTime.hour >= 12 ? "PM" : "AM";

    DateTime time = DateTime(date.year, date.month, date.day, dateTime.hour, dateTime.minute, 0);
    timeData["milliSecond"] = time.millisecondsSinceEpoch;
    timeData['object'] = time;
    return timeData;
  }

  int ipart(x) => (x > 0) ? x.floor() : x.ceil();

  double range(x) {
    double b = x / 360;
    double a = 360 * (b - ipart(b));
    if (a < 0) {
      a = a + 360;
    }
    return a;
  }

  double frac(x) {
    double a = x - x.floor();
    if (a < 0) a += 1;
    return a;
  }

  double lmst(mjd, glong) {
    double d = mjd - 51544.5;
    double t = d / 36525.0;
    double lst = range(280.46061837 + 360.98564736629 * d + 0.000387933 * t * t - t * t * t / 38710000);
    return (lst / 15.0 + glong / 15);
  }

  List<double> quad(ym, yz, yp) {
    ////  finds the parabola throuh the three points (-1,ym), (0,yz), (1, yp)//  and returns the coordinates of the max/min (if any) xe, ye//  the values of x where the parabola crosses zero (roots of the quadratic)//  and the number of roots (0, 1 or 2) within the interval [-1, 1]////  well, this routine is producing sensible answers////  results passed as array [nz, z1, z2, xe, ye]//    var nz, a, b, c, dis, dx, xe, ye, z1, z2, nz;
    List<double> quadout = [];
    double? z1, z2;
    int nz = 0;
    double a = 0.5 * (ym + yp) - yz;
    double b = 0.5 * (yp - ym);
    double c = yz;
    double xe = -b / (2 * a);
    double ye = (a * xe + b) * xe + c;
    double dis = b * b - 4.0 * a * c;
    if (dis > 0) {
      double dx = 0.5 * sqrt(dis) / a.abs();
      z1 = xe - dx;
      z2 = xe + dx;
      if (z1.abs() <= 1.0) nz += 1;
      if (z2.abs() <= 1.0) nz += 1;
      if (z1 < -1.0) z1 = z2;
    }

    quadout.add(nz.toDouble());
    quadout.add(z1!.toDouble());
    quadout.add(z2!.toDouble());
    quadout.add(xe.toDouble());
    quadout.add(ye.toDouble());
    return quadout;
  }

  dynamic hrsmin(hours) {
    ////  takes decimal hours and returns a string in hhmm format//    var hrs, h, m, dum;
    double hrs = (hours * 60 + 0.5).floor() / 60.0;
    int h = hrs.floor();
    int m = (60 * (hrs - h) + 0.5).floor();
    int dum = h * 100 + m;
    String strDum = dum.toString();
    // the jiggery pokery below is to make sure that two minutes past midnight
    // comes out as 0002 not 2. Javascript does not appear to have 'format codes'    // like C
    if (dum < 1000) strDum = "0$strDum";
    if (dum < 100) strDum = "0$strDum";
    if (dum < 10) strDum = "0$strDum";

    return strDum;
  }

  double sin_alt(iobj, mjd0, hour, glong, cglat, sglat) {
    var mjd, t, ra, dec, tau, salt;
    double rads = 0.0174532925;
    List<dynamic> objpos = [];
    mjd = mjd0 + hour / 24.0;
    t = (mjd - 51544.5) / 36525.0;
    if (iobj == 1) {
      objpos = minimoon(t);
    } else {
      objpos = minisun(t);
    }
    ra = objpos[1];
    dec = objpos[0];
    // ra = objpos[2];
    // dec = objpos[1];
    tau = 15.0 * (lmst(mjd, glong) - ra);
    salt = sglat * sin(rads * dec) + cglat * cos(rads * dec) * cos(rads * tau);
    return salt;
  }

  List<dynamic> minimoon(t) {
    double p2 = 6.283185307, arc = 206264.8062, coseps = 0.91748, sineps = 0.39778;
    double L0, L, LS, F, D, H, S, N, DL, CB, L_moon, B_moon, V, W, X, Y, Z, RHO;
    List<dynamic> mooneq = [];

    L0 = frac(0.606433 + 1336.855225 * t); // mean longitude of moon
    L = p2 * frac(0.374897 + 1325.552410 * t); //mean anomaly of Moon
    LS = p2 * frac(0.993133 + 99.997361 * t); //mean anomaly of Sun
    D = p2 * frac(0.827361 + 1236.853086 * t); //difference in longitude of moon and sun
    F = p2 * frac(0.259086 + 1342.227825 * t); //mean argument of latitude

    // corrections to mean longitude in arcsec
    DL = 22640 * sin(L);
    DL += -4586 * sin(L - 2 * D);
    DL += 2370 * sin(2 * D);
    DL += 769 * sin(2 * L);
    DL += -668 * sin(LS);
    DL += -412 * sin(2 * F);
    DL += -212 * sin(2 * L - 2 * D);
    DL += -206 * sin(L + LS - 2 * D);
    DL += 192 * sin(L + 2 * D);
    DL += -165 * sin(LS - 2 * D);
    DL += -125 * sin(D);
    DL += -110 * sin(L + LS);
    DL += 148 * sin(L - LS);
    DL += -55 * sin(2 * F - 2 * D);

    // simplified form of the latitude terms
    S = F + (DL + 412 * sin(2 * F) + 541 * sin(LS)) / arc;
    H = F - 2 * D;
    N = -526 * sin(H);
    N += 44 * sin(L + H);
    N += -31 * sin(-L + H);
    N += -23 * sin(LS + H);
    N += 11 * sin(-LS + H);
    N += -25 * sin(-2 * L + F);
    N += 21 * sin(-L + F);

    // ecliptic long and lat of Moon in rads
    L_moon = p2 * frac(L0 + DL / 1296000);
    B_moon = (18520.0 * sin(S) + N) / arc;

    // equatorial coord conversion - note fixed obliquity
    CB = cos(B_moon);
    X = CB * cos(L_moon);
    V = CB * sin(L_moon);
    W = sin(B_moon);
    Y = coseps * V - sineps * W;
    Z = sineps * V + coseps * W;
    RHO = sqrt(1.0 - Z * Z);
    double dec = (360.0 / p2) * atan(Z / RHO);
    double ra = (48.0 / p2) * atan(Y / (X + RHO));
    if (ra < 0) ra += 24;

    mooneq.add(dec);
    mooneq.add(ra);
    // mooneq[0] = dec;
    // mooneq[1] = ra;

    return mooneq;
  }

  minisun(t) {
    var p2 = 6.283185307, coseps = 0.91748, sineps = 0.39778;
    var L, M, DL, SL, X, Y, Z, RHO, ra, dec;
    List<double> suneq = [];

    M = p2 * frac(0.993133 + 99.997361 * t);
    DL = 6893.0 * sin(M) + 72.0 * sin(2 * M);
    L = p2 * frac(0.7859453 + M / p2 + (6191.2 * t + DL) / 1296000);
    SL = sin(L);
    X = cos(L);
    Y = coseps * SL;
    Z = sineps * SL;
    RHO = sqrt(1 - Z * Z);
    dec = (360.0 / p2) * atan(Z / RHO);
    ra = (48.0 / p2) * atan(Y / (X + RHO));
    if (ra < 0) ra += 24;
    suneq.add(dec);
    suneq.add(ra);
    return suneq;
  }

  double mjd() {
    double a, b;
    int month = date.month;
    int year = date.year;
    int hour = date.hour;
    int day = date.day;

    if (month <= 2) {
      month = month + 12;
      year = year - 1;
    }
    a = 10000.0 * year + 100.0 * month + day;
    if (a <= 15821004.1) {
      b = (-2 * ((year + 4716) / 4).floor() - 1179).toDouble();
    } else {
      b = ((year / 400).floor() - (year / 100).floor() + (year / 4).floor()).toDouble();
    }
    a = 365.0 * year - 679004.0;
    return a + b + (30.6001 * (month + 1)).floor() + day + hour / 24.0;
  }
}

class DeerHunting {
  DateTime date;
  double tz;
  double longitude;
  double latitude;

  DeerHunting({
    required this.date,
    required this.tz,
    required this.latitude,
    required this.longitude,
  });

  getSunMoonTime() {
    SunMoonTime? sunMoon = getSunMoonRiseSetTime(
      date: DateTime(date.year, date.month, date.day),
      latitude: latitude,
      longitude: longitude,
      tz: tz,
    );

    return sunMoon;
  }

  calculateDawnDusk(sunrise, sunset, moonrise, moonset) {
    DateTime dawn = DateTime.fromMillisecondsSinceEpoch((sunrise - (sunset - sunrise) / 2).floor());
    DateTime dusk = DateTime.fromMillisecondsSinceEpoch((sunset + (sunrise + 24 * 60 * 60 * 1000 - sunset) / 2).floor());
    DateTime moonOverhead = DateTime.fromMillisecondsSinceEpoch((moonrise + (moonset - moonrise) / 2).floor());
    DateTime moonUnderfoot = DateTime.fromMillisecondsSinceEpoch((moonset + (moonrise + 24 * 60 * 60 * 1000 - moonset) / 2).floor());

    Map<String, dynamic> deerFeedingTime = {
      "dawn": {
        "start": dawn,
        "end": DateTime.fromMillisecondsSinceEpoch(dawn.millisecondsSinceEpoch + 2 * 60 * 60 * 1000),
      },
      "dusk": {
        "start": dusk,
        "end": DateTime.fromMillisecondsSinceEpoch(dusk.millisecondsSinceEpoch - 2 * 60 * 60 * 1000),
      },
      "moonOverhead": {
        "start": DateTime.fromMillisecondsSinceEpoch(moonOverhead.millisecondsSinceEpoch - 60 * 60 * 1000),
        "end": DateTime.fromMillisecondsSinceEpoch(moonOverhead.millisecondsSinceEpoch + 60 * 60 * 1000),
      },
      "moonUnderfoot": {
        "start": DateTime.fromMillisecondsSinceEpoch(moonUnderfoot.millisecondsSinceEpoch - 60 * 60 * 1000),
        "end": DateTime.fromMillisecondsSinceEpoch(moonUnderfoot.millisecondsSinceEpoch + 60 * 60 * 1000),
      }
    };
    return deerFeedingTime;
  }

  getMajorFeedingTime(sunMoonTime) {
    Map feedDateTime = calculateDawnDusk(
      sunMoonTime.sunRiseSet()["rise"] == false ? 0 : sunMoonTime.sunRiseSet()["rise"]["milliSecond"],
      sunMoonTime.sunRiseSet()["set"] == false ? 0 : sunMoonTime.sunRiseSet()["set"]["milliSecond"],
      sunMoonTime.moonRiseSet()["rise"] == false ? 0 : sunMoonTime.moonRiseSet()["rise"]["milliSecond"],
      sunMoonTime.moonRiseSet()["set"] == false ? 0 : sunMoonTime.moonRiseSet()["set"]["milliSecond"],
    );

    DateTime moonOverHeadStart = feedDateTime["moonOverhead"]["start"];
    DateTime moonOverHeadEnd = feedDateTime["moonOverhead"]["end"];
    DateTime moonUnderfootStart = feedDateTime["moonUnderfoot"]["start"];
    DateTime moonUnderfootEnd = feedDateTime["moonUnderfoot"]["end"];

    Map overHeadStart = sunMoonTime.dateTimeObjectToHrsObject(moonOverHeadStart);
    Map overHeadEnd = sunMoonTime.dateTimeObjectToHrsObject(moonOverHeadEnd);
    Map overFootStart = sunMoonTime.dateTimeObjectToHrsObject(moonUnderfootStart);
    Map overFootEnd = sunMoonTime.dateTimeObjectToHrsObject(moonUnderfootEnd);

    var abc = overHeadStart;
    var bcd = overFootStart;
    if (overHeadStart["meridian"] == "PM") {
      abc = overHeadStart;
    }
    if (overFootStart["meridian"] == "AM") {
      bcd = overFootStart;
    }

    Map moonRTransitStart = sunMoonTime.dateTimeObjectToHrsObject(abc['object'].subtract(const Duration(hours: 1)));
    Map moonRTransitEnd = sunMoonTime.dateTimeObjectToHrsObject(abc['object'].add(const Duration(hours: 1)));

    Map moonTransitStart = sunMoonTime.dateTimeObjectToHrsObject(bcd['object']);
    Map moonTransitEnd = sunMoonTime.dateTimeObjectToHrsObject(bcd['object'].add(const Duration(hours: 2)));

    return {
      "1": {"start": moonRTransitStart, "end": moonRTransitEnd},
      "2": {"start": moonTransitStart, "end": moonTransitEnd},
    };
  }

  getMinorFeedingTime() {
    SunMoonTime sunMoonTime = getSunMoonTime();

    Map moonRiseStart = sunMoonTime.dateTimeObjectToHrsObject(sunMoonTime.moonRiseSet()["rise"]['object']);
    Map moonRiseEnd = sunMoonTime.dateTimeObjectToHrsObject(sunMoonTime.moonRiseSet()["rise"]['object'].add(const Duration(hours: 2)));

    Map moonSetStart = sunMoonTime.dateTimeObjectToHrsObject(sunMoonTime.moonRiseSet()["set"]['object']);
    Map moonSetEnd = sunMoonTime.dateTimeObjectToHrsObject(sunMoonTime.moonRiseSet()["set"]['object'].add(const Duration(hours: 2)));

    return {
      "1": {"start": moonRiseStart, "end": moonRiseEnd},
      "2": {"start": moonSetStart, "end": moonSetEnd},
    };
  }
}

class MoonWidget extends StatelessWidget {
  ///DateTime to show.
  ///Even hour, minutes, and seconds are calculated for MoonWidget
  final DateTime date;

  ///Decide the container size for the MoonWidget
  final double size;

  ///Resolution will be the moon radius.
  ///Large resolution needs more math operation makes widget heavy.
  ///Enter a small number if it is sufficient to mark it small,
  ///such as an icon or marker.
  final double resolution;

  ///Color of light side of moon
  final Color moonColor;

  ///Color of dark side of moon
  final Color earthshineColor;

  const MoonWidget({
    Key? key,
    required this.date,
    this.size = 36,
    this.resolution = 96,
    this.moonColor = Colors.amber,
    this.earthshineColor = Colors.black87,
  }) : super(key: key);

  @override
  Widget build(BuildContext context) {
    return SizedBox(
      width: size,
      height: size,
      child: Transform.scale(
        scale: size / (resolution * 2),
        child: CustomPaint(
          painter: MoonPainter(moonWidget: this),
        ),
      ),
    );
  }
}

class MoonPainter extends CustomPainter {
  MoonWidget moonWidget;
  final Paint paintDark = Paint();
  final Paint paintLight = Paint();
  final MoonPhase moon = MoonPhase();

  MoonPainter({required this.moonWidget});

  @override
  void paint(Canvas canvas, Size size) {
    double radius = moonWidget.resolution;

    int width = radius.toInt() * 2;
    int height = radius.toInt() * 2;
    double phaseAngle = moon.getPhaseAngle(moonWidget.date);

    double xcenter = 0;
    double ycenter = 0;

    try {
      paintLight.color = moonWidget.moonColor;
      //달의 색깔로 전체 원을 그린다
      canvas.drawCircle(const Offset(0, 1), radius, paintLight);
    } catch (e) {
      radius = min(width, height) * 0.4;
      paintLight.color = moonWidget.moonColor;
      Rect oval = Rect.fromLTRB(xcenter - radius, ycenter - radius, xcenter + radius, ycenter + radius);
      canvas.drawOval(oval, paintLight);
    }

    ///위상각은 태양 - 달 - 지구의 각도다.
    ///따라서 0 = full phase, 180 = new
    ///우리가 필요한 것은 일출 터미네이터의 위치 각도(태양 - 지구 - 달)다.
    ///위상각과 반대 방향이기 때문에 변환해야한다.

    double positionAngle = pi - phaseAngle;
    if (positionAngle < 0.0) {
      positionAngle += 2.0 * pi;
    }

    // myPhaseAngle = phaseAngle;
    //이제 어두운 면을 그려야 한다.
    paintDark.color = moonWidget.earthshineColor;

    double cosTerm = cos(positionAngle);

    double rsquared = radius * radius;
    double whichQuarter = ((positionAngle * 2.0 / pi) + 4) % 4;

    for (int j = 0; j < radius; ++j) {
      double rrf = sqrt(rsquared - j * j);
      double rr = rrf;
      double xx = rrf * cosTerm;
      double x1 = xcenter - (whichQuarter < 2 ? rr : xx);
      double w = rr + xx;
      canvas.drawRect(Rect.fromLTRB(x1, ycenter - j, w + x1, ycenter - j + 2), paintDark);
      canvas.drawRect(Rect.fromLTRB(x1, ycenter + j, w + x1, ycenter + j + 2), paintDark);
    }
  }

  @override
  bool shouldRepaint(CustomPainter oldDelegate) {
    return false;
  }
}

class MoonPhase {
  final deg2rad = pi / 180;

  // convert degrees to a valid angle:
  double angle(double deg) {
    while (deg >= 360.0) {
      deg -= 360.0;
    }
    while (deg < 0.0) {
      deg += 360.0;
    }
    return deg * deg2rad;
  }

  // Return the phase angle for the given date, in RADIANS.
  // Equation from Meeus eqn. 46.4.
  double getPhaseAngle(DateTime date) {
    // Time measured in Julian centuries from epoch J2000.0:
    DateTime tEpoch = DateTime(2000, 1, 1, 12);
    double t = (decimalYears(date) - decimalYears(tEpoch)) / 100.0;
    double t2 = t * t;
    double t3 = t2 * t;
    double t4 = t3 * t;

    // Mean elongation of the moon:
    double D = angle(297.8502042 + 445267.1115168 * t - 0.0016300 * t2 + t3 / 545868 + t4 / 113065000);
    // Sun's mean anomaly:
    double M = angle(357.5291092 + 35999.0502909 * t - 0.0001536 * t2 + t3 / 24490000);
    // Moon's mean anomaly:
    double mPrime = angle(134.9634114 + 477198.8676313 * t + 0.0089970 * t2 - t3 / 3536000 + t4 / 14712000);

    return (angle(180 -
        (D / deg2rad) -
        6.289 * sin(mPrime) +
        2.100 * sin(M) -
        1.274 * sin(2 * D - mPrime) -
        0.658 * sin(2 * D) -
        0.214 * sin(2 * mPrime) -
        0.110 * sin(D)));
  }

  double decimalYears(DateTime date) {
    return date.millisecondsSinceEpoch.toDouble() / 365.242191 / (24 * 60 * 60 * 1000);
  }

  double getTimeAsDecimalDay(DateTime date) {
    return date.millisecondsSinceEpoch.toDouble() / (24 * 60 * 60 * 1000);
  }
}
