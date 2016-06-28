## Raytraced Black Hole

![Screenshot](screenshot.png?raw=true "Screenshot")

Second assignment of Computer Graphics and Image Processing in spring 2015.

### Compiling and running

* Run with **run.sh**

  or

* Compile with g++
```
g++ $_FILENAME -w -g -W -Wall -o output -lGL -lglut -lm
```

### Original task
Második feladat

Készítsen sugárkövető programot, amely 100 mm oldalhosszúságú, egyik oldalán
nyitott dobozban lévő durván tesszellált, optikailag sima, tükröző arany
tóruszt, és egy Schwarzschild típusú (azaz töltés nélküli és nem forgó) Föld
tömegű fekete lyukat tartalmaz. A doboz oldalai diffúz+spekuláris
visszaverődésűek, a diffúz tag procedúrálisan textúrázott. A színteret
ambiensfény és legalább egy pontfényforrás világítja meg.
Élhet azzal az egyszerűsítéssel, hogy a fekete lyuk a fényt csak az első visszaverődés
után görbíti, azaz a pont fényforrás és a megvilágított felület között nem,
de a rücskös és sima felületek között, valamint a felületek és a szem között
igen (a másik alternatíva a fényforrásból kilépő fény görbítésére a fotontérkép
módszer).
A fény elgörbülését kellő tudásszomj esetén Einstein téregyenletének
Schwarzschildtól származó megoldásával, egyébként pedig az ekvivalenciaelv
és a klasszikus newtoni gravitációs formula felhasználásával, sugármasírozással
kell szimulálni. A rendelkezésre álló CPU idő 120 sec, szükség esetén kisebb
bontású kép átmintavételezésével, távolságtól függő lépés-nagysággal és
befoglaló térfogatokkal lehet a programot gyorsítani.

A szükséges fizikai állandók: A fény sebessége: 300 000 km/sec; A föld
kerülete: 40 000 km; A nehézségi gyorsulás a föld felszínén: 10 m/sec^2.

Az arany törésmutatója és kioltási tényezője:
r: 0.17/3.1, g: 0.35/2.7, b: 1.5/1.9
A falak diffúz és spekuláris visszaverődési tényezője és shininess paramétere,
valamint a textúra szabadon megválasztható.

Beadási határidő: 2015. 04. 12. 23:59
