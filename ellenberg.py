import csv
import statistics as stat

'''Dieses Programm ordnet den Spezies ihre Zeigerwerte zu. 
    Dazu muss es in dem selben Ordner liegen wie die aufgenommen Spezies ('daten.CSV'), 
    die ellenbergtabelle ('ellenb.CSV')und ein ergebnisfile ('zeigerwerte.CSV') liegen. 
    Das Ergebnisfile wird jedes mal wenn das Programm läuft neu überschrieben.

    Die aufgenommenen Spezies müssen in einer Spalte in daten.CSV gespeichert sein. 
    Wenn in zeigerwerte.CSV 'kein Eintrag' steht, kann das zwei Gründe haben. 
    Entweder gibt es tatsächlich keinen Tabellen eintrag für diese Spezies, oder 
    es gibt einen Tabellen eintrag unter einem anderen (alten) Namen.
    Sollte das der Fall sein können der alte und der neue Name in newoldnames
    eingetragen werden ('Name aus daten.CSV': 'Name aus ellenb.CSV'). 
    Wenn in zeigerwerte.CSV 'zuviele Einträge' steht liegt es daran, dass es 
    mehrere Einträge für subspezies gibt. Florian hat gesagt dieser Fall sei so zu
    behandeln wie 'kein Eintrag'.
    '''
class Ellenberg:
    def __init__(self):
        self.newoldnames = {'Scorzoneroides autumnalis': 'Leontodon autumnalis',
                    'Hypochaeris maculata L.':  'Hypochoeris maculata',
                    'Hypochaeris radicata': 'Hypochoeris radicata',
                    'Carex arenaria agg.': 'Carex arenaria',
                    'Achillea millefolium agg.':'Achillea millefolium',
                    'Pimpinella saxifraga agg.':'Pimpinella saxifraga',
                    'Galium mollugo agg.':'Galium mollugo',
                    'Anthoxanthum odoratum agg.':'Anthoxanthum odoratum'}

        self.headerlist = ['Spezies', 'Licht', 'Temperatur', 'Kontinentalitaet', 'Feuchte', 'Reaktion', 'Stickstoff']

    def make_speciesdict(self):
        '''Ließt die zeigerwerte aus ellenb.CSV und speichert sie in als Liste
        die einer Spezies zugeordnet ist in einem dictionary'''
        speciesdict = {}
        with open('ellenb.CSV', newline='') as csvfile:
            ellenberg = csv.reader(csvfile, delimiter=';')
            for row in ellenberg:
                speciesdict[row[1]]= [row[2],row[3],row[4],row[5],row[6],row[7]]

        return speciesdict


    def change_names_to_old(self, data):
        '''Ändert die Einträge in einer Liste entsprechend des dictionarys newoldnames
            wenn ein Eintrag als keyword in newoldnames vorkommt wird er zu dem Argument geändert'''
        change = self.newoldnames
        for i in change:
            for j in range(len(data)):
                if i ==data[j]:
                    data[j] = change[i]
                    #print('i: '+data[j])


    def change_names_to_new(self, data):
        '''Ändert die Einträge in einer Liste entsprechend des dictionarys newoldnames
            wenn ein Eintrag als argument in newoldnames vorkommt wird er zu dem Keyword geändert'''
        change = self.newoldnames
        for i in change:
            for j in range(len(data)):
                if change[i]== data[j][0]:
                    data[j][0] = i
           

    def write_results(self, results, path):
        '''Schreibt die liste in der sich die Ergebnisse befinden in zeigervalues.CSV'''
        with open (path,'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=';')
            writer.writerow(self.headerlist)
            writer.writerows(results)

    def append_to_list(self, toappend1,toappend2, result):
        '''toappend1 ist ein string, toappend2 ist eine liste von strings.
        toappend1 und alle einträge in toappend2 werden in einer neuen liste (result)
        gespeichert'''
        result.append(toappend1)
        for i in toappend2:
            result.append(i)
        return result


    def get_zeigervalues(self, speciesdict, daten):
        '''iteriert durch die liste mit den Speziesnamen (daten) und
        überprüft ob es dazu einen Eintrag in der ellenbergtabelle (speciesdict) gibt.
        Wenn es einen Eintrag gibt werden die Daten in einer Liste hinter dem Namen gespeichert.
        Wenn es keine Daten gibt wird zwischen dem Fall 'kein Eintrag' und 'zu viele Einträge' unterschieden'''
        results = []
        for i in daten:
            data = []
            smalldict = {}
            for j in speciesdict:
                if i == j:
                    self.append_to_list(i,speciesdict[j],data)
                    break
                elif i in j:
                    smalldict[j]= speciesdict[j]

            if len(smalldict) > 0:
                l = ['zuviele Eintraege','' ,'' ,'' ,'' ,'']
                self.append_to_list(i,l,data)
            if len(data)==0:
                l = ['kein Eintrag','' ,'' ,'' ,'' ,'']
                self.append_to_list(i,l,data)
            results.append(data)


        return results


    """ Diese Funktion ändern, dass ich nicht immer ein dummes csv file brauche, sondern einfach direkt mein array r
        reinladen kann """
    def make_daten(self):
        daten= []
        with open ('daten.CSV', newline='') as csvfile1:
            data = csv.reader(csvfile1, delimiter = ';')
            for row in data:
                spec = row[0]
                daten.append(spec)
        return daten
    

def do_all(column_array_species, path):
    new_ellenberg = Ellenberg()
    
    """ Ersetze das Array Daten zu meinem, was ich reinlade
        #daten = new_ellenberg.make_daten()
    """
    speciesdict = new_ellenberg.make_speciesdict()
    new_ellenberg.change_names_to_old(column_array_species)
    results = new_ellenberg.get_zeigervalues(speciesdict, column_array_species)
    new_ellenberg.change_names_to_new(results)

    licht = {"val": 0, "count": 0, "mean": None, "dev": 0}
    tempt = {"val": 0, "count": 0, "mean": None, "dev": 0}
    kont  ={"val": 0, "count": 0, "mean": None, "dev": 0}
    feuchte = {"val": 0, "count": 0, "mean": None, "dev": 0}
    reaktion = {"val": 0, "count": 0, "mean": None, "dev": 0}
    stick = {"val": 0, "count": 0, "mean": None, "dev": 0}
    
    
    array = [licht, tempt, kont, feuchte, reaktion, stick]
    
    for element in results:

        if element[1] != 'kein Eintrag':
            if element[1] != 'zuviele Eintraege':
                for counter, value in enumerate(element[1:]):
                    if value != "x" and value != "":
                        array[counter]["val"] += int(value)
                        array[counter]["count"] += 1
    
    for element in array:
        element["mean"] = element["val"]/element["count"]
        
        
    for element in results:
        if element[1] != 'kein Eintrag':
            if element[1] != 'zuviele Eintraege':
                for counter, value in enumerate(element[1:]):
                    if value != "x" and value != "":
                        array[counter]["dev"] += (int(value)-array[counter]["mean"])**2
                        
    for element in array:
        element["dev"] = (element["dev"]/element["count"])**(0.5)
                        
    new_ellenberg.write_results(results, path)
    
    return array












    