import 'dart:math';

import 'package:flutter/cupertino.dart';

class SunMoonDirection {
  String dateTimeString;
  double longitude;
  double latitude;
  DateTime dateTime = DateTime.now();
  Map<String, Map<String, dynamic>> Positions = {};
  SunMoonDirection({
    required this.latitude,
    required this.longitude,
    required this.dateTimeString,
  }) {
    dateTime = DateTime.parse(dateTimeString);
  }

  double Div(double a, double b) {
    return (a - a % b) / b;
  }

  double Rev(double number) {
    return number - (number / 360.0).floor() * 360;
  }

  double Rev2(number) {
    double x = number - (number / 360.0).floor() * 360;
    if (x > 180) x = x - 360;
    return x;
  }

  double toRad(b) {
    return b * pi / 180;
  }

  double toDeg(b) {
    return 180 * b / pi;
  }

  double daynumber(int day, int month, int year, int hour, int minute, int second) {
    double d = 367 * year - Div((7 * (year + (Div((month + 9), 12)))), 4) + Div((275 * month.toDouble()), 9) + day - 730530;
    d = d + hour / 24 + minute / (60 * 24) + second / (24 * 60 * 60);
    return d;
  }

  double ElevationRefraction(double elGeometric) {
    double El_observed;
    double x, a0, a1, a2, a3, a4;
    a0 = 0.58804392;
    a1 = -0.17941557;
    a2 = 0.29906946E-1;
    a3 = -0.25187400E-2;
    a4 = 0.82622101E-4;

    El_observed = elGeometric;

    x = (elGeometric + 0.589).abs();
    double refraction = (a0 + a1 * x + a2 * x * x + a3 * x * x * x + a4 * x * x * x * x).abs();

    if (elGeometric > 10.2) {
      El_observed = elGeometric + 0.01617 * (cos(toRad(elGeometric.abs())) / sin(toRad(elGeometric.abs())));
    } else {
      El_observed = elGeometric + refraction;
    }
    return (El_observed);
  }

  double Elevation2(double SatLon, double SiteLat, double SiteLon, double Height_over_ocean) {
    double Rstation, f, r_eq, r_sat, Ra, Rz, alfa_rx, alfa_ry, refraction, alfa_rz, alfa_r_north, alfa_r_zenith, El_geometric, El_observed;
    double x, a0, a1, a2, a3, a4;
    a0 = 0.58804392;
    a1 = -0.17941557;
    a2 = 0.29906946E-1;
    a3 = -0.25187400E-2;
    a4 = 0.82622101E-4;

    f = (1 / 298.257);
    r_sat = 42164.57;
    r_eq = 6378.14;

    Rstation = r_eq / (sqrt(1 - f * (2 - f) * sin(toRad(SiteLat)) * sin(toRad(SiteLat))));
    Ra = (Rstation + Height_over_ocean) * cos(toRad(SiteLat));
    Rz = Rstation * (1 - f) * (1 - f) * sin(toRad(SiteLat));
    // alfa_r = r_sat - Rstation;

    alfa_rx = r_sat * cos(toRad(SatLon - SiteLon)) - Ra;
    alfa_ry = r_sat * sin(toRad(SatLon - SiteLon));
    alfa_rz = -Rz;

    alfa_r_north = -alfa_rx * sin(toRad(SiteLat)) + alfa_rz * cos(toRad(SiteLat));
    alfa_r_zenith = alfa_rx * cos(toRad(SiteLat)) + alfa_rz * sin(toRad(SiteLat));

    El_geometric = toDeg(atan2(alfa_r_zenith, sqrt(alfa_r_north * alfa_r_north + alfa_ry * alfa_ry)));

    x = (El_geometric + 0.589).abs();
    refraction = (a0 + a1 * x + a2 * x * x + a3 * x * x * x + a4 * x * x * x * x).abs();

    if (El_geometric > 10.2) {
      El_observed = El_geometric + 0.01617 * (cos(toRad(El_geometric.abs())) / sin(toRad(El_geometric.abs())));
    } else {
      El_observed = El_geometric + refraction;
    }
    return (El_observed);
  }

  dynamic formatnumber(double num, places) {
    double i, integer;
    List ss = [];
    double numcopy = num;
    var a = pow(10, places.abs() == places ? places : 2);
    String strOP = ((num * a).round() / a).toString();
    if (num < 0) strOP = " $strOP";
    ss = strOP.split(".");
    if ((num != 0) && (ss.length > 1)) {
      integer = double.parse(ss[0].trim());
      double decimals = double.parse(ss[1].trim());
      if (decimals.toString().length < places) {
        String addzeroes = "0";
        for (i = 0; i < (places - decimals.toString().length - 1); i++) {
          addzeroes += "0";
        }
        strOP = "$integer.$decimals$addzeroes";
      }
    } else {
      String decimals = "0";
      for (i = 0; i < places - 1; i++) {
        decimals += "0";
      }
      strOP = "$numcopy.$decimals";
    }
    if (numcopy >= 0) strOP = " $strOP";
    strOP = strOP.trim();
    return strOP;
  }

  String formatvalue(input, int rsize) {
    String invalid = "**************************";
    String nines = "999999999999999999999999";
    String strin = "$input";
    double fltin = double.parse(strin);
    if (strin.length <= rsize) return strin;
    if (strin.contains("e") || fltin > double.parse("${nines.substring(0, rsize)}.4")) {
      return invalid.substring(0, rsize);
    }
    var rounded = "$fltin + ${(fltin - double.parse(strin.substring(0, rsize)))}";
    return rounded.substring(0, rsize);
  }

  List getsatelliteposition(double InputAzimuth, double latitude, double longitude) {
    String? SatteliteLon, SatLon;
    double X, Y, Z, R, S, p, Azimuth;
    R = 6378.14; //' km
    S = 42164.57; //' km

    double Lat = 1 * latitude;
    double Lon = 1 * longitude;

    Azimuth = InputAzimuth;
    Z = (Lon * pi / 180) - atan(tan(Azimuth * pi / 180) * sin(Lat * pi / 180));
    Z = Rev(Z * 180 / pi);
    if ((Azimuth > 90) && (Azimuth < 270)) {
      if (Z < 0) {
        Z = 360 - Z;
        SatteliteLon = formatnumber(Z, 5) + " °W*";
        SatLon = Z.toString();
      }
      if (Z > 0) {
        if (Z > 180) {
          Z = 360 - Z;
          SatteliteLon = formatnumber(Z, 5) + " °W (" + formatnumber((360 - Z), 3) + " °E)";
          SatLon = (-Z).toString();
        } else {
          SatteliteLon = formatnumber(Z, 5) + " °E";
          SatLon = Z.toString();
        }
      }
    } else {
      SatteliteLon = "Not Available";
      SatLon = Z.toString();
    }

    return [SatteliteLon, SatLon];
  }

  List sun_angles(double d, double SiteLon, double SiteLat) {
    double HourAngle, SIDEREALTIME;
    double w, a, e, M, L, oblecl, E, x, y, r, v, sunlon, z, xequat, yequat, zequat, RA, Decl, GMST0, UT;
    double SunElevation, xhor, yhor, zhor, GeometricElevation;
    List angles = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    w = 282.9404 + 4.70935E-5 * d;
    a = 1;
    e = 0.016709 - 1.151E-9 * d;
    M = 356.0470 + 0.9856002585 * d;
    oblecl = 23.4393 - 3.563E-7 * d;
    L = w + Rev(M);
    L = Rev(L);
    E = M + (180 / pi) * e * sin(toRad(M)) * (1 + e * cos(toRad(M)));
    E = Rev(E);
    x = a * cos(toRad(E)) - e;
    y = a * sin(toRad(Rev(E))) * sqrt(1 - e * e);
    r = sqrt(x * x + y * y);
    v = toDeg(atan2(y, x));
    sunlon = Rev(v + w);

    x = r * cos(toRad(sunlon));
    y = r * sin(toRad(sunlon));
    z = 0;

    xequat = x;
    yequat = y * cos(toRad(oblecl)) + z * sin(toRad(oblecl));
    zequat = y * sin(toRad(oblecl)) + z * cos(toRad(oblecl));

    RA = Rev(toDeg(atan2(yequat, xequat)));
    Decl = toDeg(atan2(zequat, sqrt(xequat * xequat + yequat * yequat)));

    GMST0 = (L + 180);

    UT = d - (d).floor();
    SIDEREALTIME = GMST0 + UT * 360 + SiteLon;
    HourAngle = SIDEREALTIME - RA;

    x = cos(HourAngle * pi / 180) * cos(Decl * pi / 180);
    y = sin(HourAngle * pi / 180) * cos(Decl * pi / 180);
    z = sin(Decl * pi / 180);

    xhor = x * sin(SiteLat * pi / 180) - z * cos(SiteLat * pi / 180);
    yhor = y;
    zhor = x * cos(SiteLat * pi / 180) + z * sin(SiteLat * pi / 180);

    double EarthLat = toDeg(atan2(z, sqrt(x * x + y * y)));
    double EarthLon = Rev((0 * 180 + RA - GMST0 - UT * 360).toDouble());

    SunElevation = toDeg(asin(zhor));
    GeometricElevation = SunElevation;
    SunElevation = ElevationRefraction(SunElevation);
    double SunAzimuth = toDeg(atan2(yhor, xhor));

    angles[0] = SunElevation;
    // angles.insert(0, SunElevation);
    angles[1] = SunAzimuth + 180;
    // angles.insert(1, SunAzimuth + 180);
    angles[2] = Decl;
    // angles.insert(2, Decl);
    angles[3] = sunlon;
    // angles.insert(3, sunlon);
    angles[4] = RA;
    // angles.insert(4, RA);
    angles[5] = Rev(GMST0);
    // angles.insert(5, Rev(GMST0));
    angles[6] = Rev(M);
    // angles.insert(6, Rev(M));
    angles[7] = Rev(w);
    // angles.insert(7, Rev(w));
    angles[8] = Rev(e);
    // angles.insert(8, Rev(e));
    angles[9] = Rev(oblecl);
    // angles.insert(9, Rev(oblecl));
    angles[10] = GeometricElevation;
    // angles.insert(10, GeometricElevation);
    angles[11] = L;
    // angles.insert(11, L);
    // angles.insert(12, 0);
    // angles.insert(13, 0);

    if (SiteLat < 0) {
      angles[14] = Rev((360 - HourAngle).toDouble());
      // angles.insert(14, Rev((360 - HourAngle).toDouble()));
    } else {
      angles[14] = Rev((HourAngle - 180).toDouble());
      // angles.insert(14, Rev((HourAngle - 180).toDouble()));
    }
    // angles.insert(15, 0);
    // angles.insert(16, 0);
    // angles.insert(17, 0);
    // angles.insert(18, 0);
    angles[19] = Rev(EarthLon);
    // angles.insert(19, Rev(EarthLon));
    angles[20] = EarthLat;
    // angles.insert(20, EarthLat);
    angles[21] = Rev2(HourAngle);
    // angles.insert(21, Rev2(HourAngle));

    return (angles);
  }

  List moon_angles(double d, double SiteLon, double SiteLat) {
    double HourAngle, SIDEREALTIME;
    double w, a, e, M, L, N, oblecl, E, x, y, r, v, sunlon, z, xequat, yequat, zequat, RA, Decl, GMST0, UT, SIDTIME, HA;
    double xhor, yhor, zhor, GeometricElevation;
    double E0, E1, xeclip, yeclip, zeclip, Lm, Ls, Ms, Mm, D, F;
    double P_lon1, P_lon2, P_lon3, P_lon4, P_lon5, P_lon6, P_lon7, P_lon8, P_lon9, P_lon10, P_lon11, P_lon12;
    double P_lat1, P_lat2, P_lat3, P_lat4, P_lat5, P_lat, P_lon, P_moondistance, Moon_RA, Moon_Decl;
    double xh, yh, zh, MoonAzimuth, MoonElevation, Iterations, E_error, Ebeforeit, Eafterit, E_ErrorBefore;
    List angles = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    List sunangles = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

    N = 125.1228 - 0.0529538083 * d;
    double i = 5.1454;

    w = 318.0634 + 0.1643573223 * d;
    a = 60.2666;
    e = 0.054900;
    M = 115.3654 + 13.0649929509 * d;

    w = Rev(w);
    M = Rev(M);
    N = Rev(N);

    E = M + (180 / pi) * e * sin(toRad(M)) * (1 + e * cos(toRad(M)));
    E = Rev(E);

    Ebeforeit = E;

    Iterations = 0;
    E_error = 9;

    while ((E_error > 0.0005) && (Iterations < 20)) {
      Iterations = Iterations + 1;
      E0 = E;
      E1 = E0 - (E0 - (180 / pi) * e * sin(toRad(E0)) - M) / (1 - e * cos(toRad(E0)));
      E = Rev(E1);

      Eafterit = E;

      if (E < E0) {
        E_error = E0 - E;
      } else {
        E_error = E - E0;
      }

      if (E < Ebeforeit) {
        E_ErrorBefore = Ebeforeit - E;
      } else {
        E_ErrorBefore = E - Ebeforeit;
      }

      // window.status = "Iterations=" + Iterations + " Moon eccentric anomaly error before iterations=" + formatvalue(E_ErrorBefore, 7) + "_deg after iterations=" + E_error + "_deg";
      // if (Iterations > 10) window.status = "Number of Iterations more than 10=" + Iterations + " Ebefore=" + Ebeforeit + " Eafter=" + Eafterit + " E_errorbefore=" + formatvalue(E_ErrorBefore, 7) + " E_errorafter" + E_error;
    }

    x = a * (cos(toRad(E)) - e);
    y = a * sin(toRad(Rev(E))) * sqrt(1 - e * e);
    r = sqrt(x * x + y * y);
    v = toDeg(atan2(y, x));

    sunlon = Rev(v + w);

    x = r * cos(toRad(sunlon));
    y = r * sin(toRad(sunlon));
    z = 0;

    xeclip = r * (cos(toRad(N)) * cos(toRad(v + w)) - sin(toRad(N)) * sin(toRad(v + w)) * cos(toRad(i)));
    yeclip = r * (sin(toRad(N)) * cos(toRad(v + w)) + cos(toRad(N)) * sin(toRad(v + w)) * cos(toRad(i)));
    zeclip = r * sin(toRad(v + w)) * sin(toRad(i));

    double moon_longitude = Rev(toDeg(atan2(yeclip, xeclip)));
    double moon_latitude = toDeg(atan2(zeclip, sqrt(xeclip * xeclip + yeclip * yeclip)));

    sunangles = sun_angles(d, SiteLon, SiteLat);
    Ls = sunangles[11]; // Suns mean longitude er feil
    Ms = sunangles[6]; // Suns mean anomaly
    Mm = Rev(M); // Moons mean anomaly
    Lm = Rev(N + w + M); // moon mean longitude
    D = Rev(Lm - Ls); //Moons mean elongation
    F = Rev(Lm - N); //Moons argument of latitude

    P_lon1 = -1.274 * sin(toRad(Mm - 2 * D)); //  (Evection)
    P_lon2 = 0.658 * sin(toRad(2 * D)); //    (Variation)
    P_lon3 = -0.186 * sin(toRad(Ms)); //    (Yearly equation)
    P_lon4 = -0.059 * sin(toRad(2 * Mm - 2 * D));
    P_lon5 = -0.057 * sin(toRad(Mm - 2 * D + Ms));
    P_lon6 = 0.053 * sin(toRad(Mm + 2 * D));
    P_lon7 = 0.046 * sin(toRad(2 * D - Ms));
    P_lon8 = 0.041 * sin(toRad(Mm - Ms));
    P_lon9 = -0.035 * sin(toRad(D)); //      (Parallactic equation)
    P_lon10 = -0.031 * sin(toRad(Mm + Ms));
    P_lon11 = -0.015 * sin(toRad(2 * F - 2 * D));
    P_lon12 = 0.011 * sin(toRad(Mm - 4 * D));

    P_lat1 = -0.173 * sin(toRad(F - 2 * D));
    P_lat2 = -0.055 * sin(toRad(Mm - F - 2 * D));
    P_lat3 = -0.046 * sin(toRad(Mm + F - 2 * D));
    P_lat4 = 0.033 * sin(toRad(F + 2 * D));
    P_lat5 = 0.017 * sin(toRad(2 * Mm + F));

    P_lon = P_lon1 + P_lon2 + P_lon3 + P_lon4 + P_lon5 + P_lon6 + P_lon7 + P_lon8 + P_lon9 + P_lon10 + P_lon11 + P_lon12;
    P_lat = P_lat1 + P_lat2 + P_lat3 + P_lat4 + P_lat5;
    P_moondistance = -0.58 * cos(toRad(Mm - 2 * D)) - 0.46 * cos(toRad(2 * D));

    moon_longitude = moon_longitude + P_lon;
    moon_latitude = moon_latitude + P_lat;
    r = r + P_moondistance;

    xh = r * cos(toRad(moon_longitude)) * cos(toRad(moon_latitude));
    yh = r * sin(toRad(moon_longitude)) * cos(toRad(moon_latitude));
    zh = r * sin(toRad(moon_latitude));

    xequat = xh;
    yequat = yh * cos(toRad(sunangles[9])) - zh * sin(toRad(sunangles[9]));
    zequat = yh * sin(toRad(sunangles[9])) + zh * cos(toRad(sunangles[9]));
    Moon_RA = Rev(toDeg(atan2(yh, xh)));
    Moon_Decl = toDeg(atan2(zh, sqrt(xh * xh + yh * yh)));

    Moon_RA = Rev(toDeg(atan2(yequat, xequat)));
    Moon_Decl = toDeg(atan2(zequat, sqrt(xequat * xequat + yequat * yequat)));

    GMST0 = (Ls + 180);
    UT = d - (d).floor();

    SIDEREALTIME = GMST0 + UT * 360 + SiteLon;
    HourAngle = SIDEREALTIME - Moon_RA;

    x = cos(HourAngle * pi / 180) * cos(Moon_Decl * pi / 180);
    y = sin(HourAngle * pi / 180) * cos(Moon_Decl * pi / 180);
    z = sin(Moon_Decl * pi / 180);

    xhor = x * sin(SiteLat * pi / 180) - z * cos(SiteLat * pi / 180);
    yhor = y;
    zhor = x * cos(SiteLat * pi / 180) + z * sin(SiteLat * pi / 180);

    MoonElevation = toDeg(asin(zhor));

    MoonElevation = MoonElevation - toDeg(asin(1 / r * cos(toRad(MoonElevation))));

    GeometricElevation = MoonElevation;
    MoonElevation = ElevationRefraction(MoonElevation);
    MoonAzimuth = toDeg(atan2(yhor, xhor));
    angles[0] = MoonElevation;

    angles[1] = MoonAzimuth + 180;

    double mpar = asin(1 / r);
    double gclat, topRA, gcHA, topDecl, g, rho;

    double EarthLon = toDeg(atan2(y, x));
    double EarthLat = toDeg(atan2(z, sqrt(x * x + y * y)));

    EarthLon = toDeg(atan2(yhor, xhor));
    gclat = (SiteLat * (pi / 180) - 0.1924 * (pi / 180) * sin(2 * SiteLat * pi / 180));
    rho = 0.99833 + 0.00167 * cos(2 * SiteLat * pi / 180);

    g = atan(tan(gclat) / cos(HourAngle * pi / 180));

    topRA = Moon_RA * pi / 180 - (mpar * rho * cos(gclat) * sin(HourAngle * pi / 180) / cos(Moon_Decl * pi / 180));
    topDecl = Moon_Decl * pi / 180 - (mpar * rho * sin(gclat) * sin(g - Moon_Decl * (pi / 180)) / sin(g));

    gcHA = (SIDEREALTIME * (pi / 180) - topRA);

    EarthLon = Rev((0 * 180 + topRA * 180 / pi - GMST0 - UT * 360).toDouble());

    angles[2] = Moon_Decl;
    angles[3] = moon_longitude;
    angles[4] = Moon_RA;
    angles[5] = Rev(GMST0);
    angles[6] = Rev(M);
    angles[7] = Rev(w);
    angles[8] = Rev(e);
    angles[9] = Rev(0);
    angles[10] = GeometricElevation;
    angles[11] = moon_latitude;
    angles[12] = MoonElevation;
    angles[13] = r;

    if (SiteLat < 0) {
      angles[14] = Rev((360 - HourAngle).toDouble());
    } else {
      angles[14] = Rev(HourAngle - 180);
    }
    if (SiteLat < 0) {
      angles[15] = Rev((360 - ((gcHA) * 180 / pi)).toDouble());
    } else {
      angles[15] = Rev(((gcHA) * 180 / pi) - 180);
    }
    angles[16] = topRA * 180 / pi;
    angles[17] = topDecl * 180 / pi;
    angles[18] = gclat * 180 / pi;
    angles[19] = EarthLon;
    angles[20] = EarthLat;

    angles[21] = Rev2(HourAngle);
    angles[22] = Rev2(((gcHA) * 180 / pi));
    return (angles);
  }

  Map<String, dynamic> CalculateAngles() {
    int Day, Month, Year, Hour, Minute, Seconds;
    var Lat, Lon, LatDir, LonDir, SunAltitude, ElevationDifference;
    var Hourangle_direction_text, Latitude_direction_text;
    double MonSatLon, SunSatLon;
    List SunSatLonArray = [];
    Map<String, Map<String, dynamic>> Angles = {};
    // double latitude, double longitude, DateTime date
    DateTime? date = dateTime;
    Angles = {
      "moon": {
        "eclipticLongitude": null,
        "eclipticLatitude": null,
        "rightAscension": null,
        "declination": null,
        "azimuth": null,
        "distance": null,
        "azimuth": null,
        "geometricElevation": null,
        "elevation": null,
        "hourAngle": null,
        "topocentric": {
          "hourAngle": "",
          "rightAscension": "",
          "declination": "",
        },
        "earth": {
          "latitude": "",
          "longitude": "",
        }
      },
      "sun": {
        "azimuth": null,
        "declination": null,
        "longitude": null,
        "rightAscension": null,
        "GMST0": null,
        "meanAnomaly": null,
        "perihelionLongitude": null,
        "eccentricity": null,
        "obliquity": null,
        "geometricElevation": null,
        "hourAngle": null,
        'earth': {
          "latitude": "",
          'longitude': '',
        }
      },
      'satellite': {
        'sun': {
          'longitude': "",
          'elevation': null,
          'elevationDiff': null,
        },
        'moon': {
          'longitude': "",
          'elevation': null,
          'elevationDiff': null,
        }
      }
    };

    DateTime newDate = date.toUtc();
    int year = newDate.year;
    int month = newDate.month;
    int day = newDate.day;
    int hour = newDate.hour;
    int minute = newDate.minute;
    int second = newDate.second;
    debugPrint("day::$day, month::$month, year::$year, hour::$hour, minute::$minute");
    double SiteLat = 1 * latitude;
    double SiteLon = 1 * longitude;

    Year = 1 * year;
    Month = 1 * month;
    Day = 1 * day;
    Hour = 1 * hour;
    Minute = 1 * minute;
    Seconds = 1 * second;

    double d = daynumber(Day, Month, Year, Hour, Minute, Seconds);

    List output_angles = sun_angles(d, SiteLon, SiteLat);
    List output_moon_angles = moon_angles(d, SiteLon, SiteLat);

    Angles['moon']?['eclipticLongitude'] = formatnumber(output_moon_angles[3], 3);
    Angles['moon']?['eclipticLatitude'] = formatnumber(output_moon_angles[11], 3);
    Angles['moon']?['rightAscension'] = formatnumber(output_moon_angles[4], 3);
    Angles['moon']?['declination'] = formatnumber(output_moon_angles[2], 3);
    Angles['moon']?['azimuth'] = formatnumber(output_moon_angles[1], 3);
    Angles['moon']?['distance'] = formatnumber(output_moon_angles[13] * 6378.14, 3);
    Angles['moon']?['geometricElevation'] = formatnumber(output_moon_angles[10], 3);
    Angles['moon']?['elevation'] = formatnumber(output_moon_angles[12], 3);

    Hourangle_direction_text = '';
    if (output_moon_angles[21] < 0) {
      Hourangle_direction_text = ' °E';
    } else {
      Hourangle_direction_text = ' °W';
    }

    Angles['moon']?['hourAngle'] = formatnumber(output_moon_angles[21], 3) + Hourangle_direction_text + ' (' + formatnumber(output_moon_angles[14], 3) + ')';

    if (1 * output_moon_angles[20] < 0) {
      Latitude_direction_text = ' °S';
    } else {
      Latitude_direction_text = ' °N';
    }

    if (output_moon_angles[22] < 0) {
      Hourangle_direction_text = ' °E';
    } else {
      Hourangle_direction_text = ' °W';
    }

    Angles['moon']?['topocentric']?['hourAngle'] = formatnumber(output_moon_angles[22], 3) + Hourangle_direction_text + ' (' + formatnumber(output_moon_angles[15], 3) + ')';
    Angles['moon']?['topocentric']?['rightAscension'] = formatnumber((1 * output_moon_angles[16]).toDouble(), 3);
    Angles['moon']?['topocentric']?['declination'] = formatnumber((1 * output_moon_angles[17]).toDouble(), 3);
    Angles['moon']?['earth']?['latitude'] = formatnumber((1 * output_moon_angles[20]).toDouble(), 3) + Latitude_direction_text;
    Angles['moon']?['earth']?['longitude'] = formatnumber((1 * output_moon_angles[19]).toDouble(), 3) + '°E ( ' + formatnumber((360 - 1 * output_moon_angles[19]).toDouble(), 3) + '°W )';

    if (double.parse(Angles['moon']!['geometricElevation']!) > -5) {
      Angles['moon']?['elevation'] = formatvalue(output_moon_angles[12], 5);
    } else {
      Angles['moon']?['elevation'] = 'N/A';
    }

    SunAltitude = output_angles[0];

    Angles['sun']?['azimuth'] = formatvalue(output_angles[1], 8);
    Angles['sun']?['declination'] = formatnumber(output_angles[2], 3);
    Angles['sun']?['longitude'] = formatnumber(output_angles[3], 3);
    Angles['sun']?['rightAscension'] = formatnumber(output_angles[4], 3);
    Angles['sun']?['GMST0'] = formatnumber(output_angles[5], 3);
    Angles['sun']?['meanAnomaly'] = formatnumber(output_angles[6], 3);
    Angles['sun']?['perihelionLongitude'] = formatnumber(output_angles[7], 3);
    Angles['sun']?['eccentricity'] = formatnumber(output_angles[8], 3);
    Angles['sun']?['obliquity'] = formatnumber(output_angles[9], 3);
    Angles['sun']?['geometricElevation'] = formatnumber(output_angles[10], 3);

    Latitude_direction_text = '';
    Hourangle_direction_text = '';
    if (1 * output_angles[20] < 0) {
      Latitude_direction_text = ' °S';
    } else {
      Latitude_direction_text = ' °N';
    }

    if (output_angles[21] < 0) {
      Hourangle_direction_text = ' °E';
    } else {
      Hourangle_direction_text = ' °W';
    }

    Angles['sun']?['hourAngle'] = formatnumber(output_angles[21], 3) + Hourangle_direction_text + ' (' + formatnumber(output_angles[14], 3) + ')';
    Angles['sun']?['earth']?['longitude'] = formatnumber(output_angles[19], 3) + '°E ( ' + formatnumber((360 - output_angles[19]).toDouble(), 3) + '°W )';
    Angles['sun']?['earth']?['latitude'] = formatnumber(output_angles[20], 3) + Latitude_direction_text;

    if (double.parse(Angles['sun']?['geometricElevation']) > -5) {
      Angles['sun']?['elevation'] = formatnumber(SunAltitude, 5);
    } else {
      Angles['sun']?['elevation'] = 'N/A';
    }

    SunSatLonArray = getsatelliteposition(output_angles[1], SiteLat, SiteLon);
    Angles['satellite']?['sun']?['longitude'] = SunSatLonArray[0];
    SunSatLon = double.parse(SunSatLonArray[1]) * 1;

    if (Angles['satellite']?['sun']?['longitude'] == "Not Available") {
      Angles['satellite']?['sun']?['elevation'] = "N/A";
      ElevationDifference = "N/A";
    } else {
      Angles['satellite']?['elevation'] = formatnumber(Elevation2(SunSatLon, SiteLat, SiteLon, 0), 5);
      ElevationDifference =
          (Angles['sun']?['elevation'] == null ? 0.0 : double.parse(Angles['sun']?['elevation']) - (Angles["satellite"]?["sun"]?["elevation"] == null ? 0.0 : double.parse(Angles["satellite"]?["sun"]?["elevation"]))).abs();
    }

    Angles['satellite']?["sun"]?["elevationDiff"] = formatnumber(0.5, 5);

    List MoonLonArray = getsatelliteposition(output_moon_angles[1], SiteLat, SiteLon);
    Angles['satellite']?['moon']?['longitude'] = MoonLonArray[0];
    MonSatLon = double.parse(MoonLonArray[1]) * 1;

    if (Angles['satellite']?['moon']['longitude'] == "Not Available") {
      Angles['satellite']?['moon']['elevation'] = "N/A";
      ElevationDifference = "N/A";
    } else {
      Angles['satellite']?['moon']['elevation'] = formatnumber(Elevation2(MonSatLon, SiteLat, SiteLon, 0), 5);
      ElevationDifference = (Angles['moon']?['elevation'] == "N/A" ? 0.0 : double.parse(Angles['moon']?['elevation']) - double.parse(Angles['satellite']?['moon']?['elevation'])).abs();
    }

    Angles["satellite"]?['moon']?['elevationDiff'] = formatnumber(ElevationDifference == "N/A" ? 0.0 : ElevationDifference, 5);
    Positions = Angles;
    return Angles;
  }

  getHuntingIndex() {
    CalculateAngles();
    return (cos(toRad(double.parse(Positions['sun']?['geometricElevation'])).abs()) + sin(toRad(double.parse(Positions['moon']?['geometricElevation']))).abs()) / 2;
  }
}
