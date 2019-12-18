# FYS3150_Project_5
Partial Differential equations; Diffusion Equation

I dette prosjektet har vi tatt for oss et geofysisk problem, vi har sett på litosfæren utenfor norskekysten der under et fenomen i en suduksjonssone, blir mantelen blir beriket med radioaktive stoffer. Dette medfører en økning av varmeprodukjsonen i dette område. Vi ønsker å finne ut om det vil være mulig å finne bevis, spor eller tegn på at et slikt fenomen har skjedd for over et giga år siden. Derfor har vi sett på utviklingen av temperatur distrubisjonen over et giga år.

For å gjøre disse simulasjonene har vi brukt paritsielle diffrensial likninger i både (1+1) og (2+1) dimensjoner. For å løse disse partisielle differnsial likningene har vi analysert den eksplisitte metoden Forward Eulers og de implisitte metodene Backward Euler og Crank Nicolson. For å ta for oss hvem av disse som er best egnet til å finne grensebetingelsene for simulasjonen.

I dette prosjektet har vi kommet frem til at Crank Nicolson er den beste metoden når det kommer til presisjon, men når det kommer til tidsforbruk vil Forward Euler være best. Videre har vi utviklet en 2 dimensjonal algoritme som gir oss fornuftige løsninger for kjente analytiske løsningers. Under simuleringene av litosfæren med lokal radioaktiv berikelse holder en ( 2 + 1 ) dimensjonal løser på grunn av symmetri, dette gjør at denne metoden passer svært godt til å simulere denne problemstillingen med.

Resultatene vi får fra simulasjonene ser veldig fornuftige og realistiske ut, og bygger meget godt på hverandre, selv om selve simulasjonstiden ble relativt lang. Resultatene tyder på at det vil være tydelige spor eller tegn på at en slik hendelse har funnet sted.
