# Alles staat in [deze](http://htmlpreview.github.io/?https://github.com/Nikdwal/finite-fields/blob/master/notebook.html) Jupyter-notebook

#### Let op de volgende zaken:

- Elk element van een veld: getallen, primitieve elementen zijn van het type FieldElement, dus NIET int of iets dergelijk. Het kan handig zijn om hun veld op te vragen met de variabele "field"

- Probeer nooit elementen van twee verschillende velden bij elkaar op te tellen of te vermenigvuldigen. Zelfs als het ene veld een uitbreidingsveld is van een ander. Bijvoorbeeld: als je een element "1" in GF(2) hebt, en een element "1" in GF(4), kun je deze twee niet met elkaar optellen. Formeel zijn deze twee objecten van het type FieldElement met een andere waarde voor de variabele "field".

- Dit betekent dus dat je best geen elementen van een veld E probeert in te vullen in een veelterm die gedefinieerd is over een ander veld F, zelfs niet als E het veld F volledig omvat. Maak in plaats daarvan een nieuwe veelterm aan over E, met coefficienten in E.
