# Geomagnetism
Java class to calculate magnetic declination, magnetic field strength, inclination etc. for any point on the earth.

## Usage:
```java
// Create instance with only location and without altitude & date
Geomagnetism gm = new Geomagnetism(123.45678, 76.54321);
System.out.printf("Declination: %f\n", gm.getDeclination());

// Create a calendar for calculation
Calendar calendar = Calendar.getInstance();
calendar.set(2017, 6, 1); // 2017 July 01
// Update instance with altitude and date apart from location
gm.calculate(98.76543, 12.34567, 1234, calendar);
System.out.printf("Declination: %f\n", gm.getDeclination());
```
