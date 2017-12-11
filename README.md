# Geomagnetism
Java class to calculate magnetic declination, magnetic field strength, inclination etc. for any point on the earth.

## Usage:
```java
// Create instance
Geomagnetism gm = new Geomagnetism(123.45678, 76.54321);
System.out.printf("Declination: %f\n", gm.getDeclination());
// Update instance
gm.calculate(98.76543, 12.34567);
System.out.printf("Declination: %f\n", gm.getDeclination());
```
