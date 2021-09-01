Questo file README.txt è stato generato in data 01/09/2021 da MATTIA ERROI


In breve: 

	Codici di calcolo per l'analisi delle interazioni plasma-parete in presenza di un divertore a metallo liquido.

	Questi codici di calcolo sono stati implementati per la realizzazione della tesi di laurea triennale in Ingegneria Energetica presso il Politecnico di Torino.

INFORMAZIONI GENERALI

1. Informazioni sull'autore

		Nome: Mattia Erroi
                
                Studente presso il Politecnico di Torino

		Istituzione (alla data di creazione): Politecnico di Torino
		Indirizzo: Corso Duca degli Abruzzi, 24
		Email personale: mattia.erroi7@gmail.com
                Email istituzionale: s256668@studenti.polito.it

		
2. Data di realizzazione  

	dal 01/03/2021 al 01/09/2021


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: 

	Open Access under Creative Commons Attribution ShareAlike CC BY-SA 

2. Links to publications that cite or use the data: 

	N/A

3. Links to other publicly accessible locations of the data: 

	N/A

4. Links/relationships to ancillary data sets: 

	N/A

5. Was data derived from another source?

	No.

6. Recommended citation for this dataset: 

	N/A

PANORAMICA DEI CONTENUTI

-cartella "codice_no_vapour_shielding"
            
          contenuti:
          -modello_1D_parete_sn.m: script MATLAB con l'analisi termica del target del divertore per le interazioni plasma-parete in assenza di vapour shielding
          -myNewton.m: script MATLAB contentente una function con il metodo di Newton per la risoluzione dell'equazione non lineare di Bergles-Rohsenow

          -PARTE1.SLDPRT: modello SOLIDWORKS per il componente del target di CuCrZr


-cartella "codice_plasma"

          contenuti:
          -modello_plasma_sn.m: script MATLAB con l'implementazione del modello di Lengyel. Si tratta di uno script ausiliario per verificare la corretta lettura dei dati della funzione di                                 radiazione e la correttezza della risoluzione del sistema non lineare del modello di Lengyel
          -myNewton_jac.m: script MATLAB contentente una function con il metodo di Newton per la risoluzione del sistema non lineare del modello di Lengyel
          -numerical_jacobian.m: script MATLAB contentente una function per il calcolo numerico della matrice jacobiana per il metodo di Newton

          -Sn_dati.dat: file contenente i dataset per la lettura dallo script di MATLAB della funzione di radiazione dello stagno

-cartella "codice_vapour_shielding"

           contenuti:
           -codice_finale_sn.m: script MATLAB con l'analisi termica del target del divertore per le interazioni plasma-parete in presenza di vapour shielding
           -myNewton.m: script MATLAB contentente una function con il metodo di Newton per la risoluzione dell'equazione non lineare di Bergles-Rohsenow
           -myNewton_jac.m: script MATLAB contentente una function con il metodo di Newton per la risoluzione del sistema non lineare del modello di Lengyel
           -numerical_jacobian.m: script MATLAB contentente una function per il calcolo numerico della matrice jacobiana per il metodo di Newton

           -Sn_dati.dat: file contenente i dataset per la lettura dallo script di MATLAB della funzione di radiazione dello stagno

-cartella "convergenza"

           contenuti:
           -convergenza_spazio.m: script MATLAB con lo studio di convergenza spaziale del codice per l'analisi termica in assenza di vapour shielding
           -myNewton.m: script MATLAB contentente una function con il metodo di Newton per la risoluzione dell'equazione non lineare di Bergles-Rohsenow
           -convergenza_tempo.m: script MATLAB con lo studio di convergenza temporale del codice per l'analisi termica in assenza di vapour shielding

ISTRUZIONI

-modello_plasma_sn.m: avviare lo scrpt su MATLAB
-modello_plasma_sn.m: avviare lo scrpt su MATLAB
-convergenza_tempo.m: avviare lo scrpt su MATLAB
-convergenza_spazio.m: avviare lo scrpt su MATLAB
-codice_finale_sn.m: avviare lo scrpt su MATLAB. Nota: alla riga 150 si può cambiare il valore della variabile "coeff_mitigazione" 


