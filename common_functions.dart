
import 'moon_phase_functions.dart';

getSunMoonRiseSetTime({
  DateTime? date,
  required double latitude,
  required double longitude,
  required double tz,
}) {
  return SunMoonTime(
    date: date ?? DateTime.now(),
    latitude: latitude,
    longitude: longitude,
    tz: tz,
  );
}

