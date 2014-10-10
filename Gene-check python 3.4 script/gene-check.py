__author__ = 'anonym'

import argparse
import io

#################### Глобальные переменные ####################

# DEBUG_ выставляйте в True при дебаге
DEBUG_ = True

codon_start_dna = 'ATG'
codon_stop_dna  = ['TAA','TAG','TGA']

codon_start_rna = 'AUG'
codon_stop_rna  = ['UAA','UAG','UGA']

RnaToProtein =\
    {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
     'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
     'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
     'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
     'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
     'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
     'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
     'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
     'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
     'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
     'UAA': '$', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
     'UAG': '$', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
     'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
     'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
     'UGA': '$', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
     'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}

ProteinToRna = {'A': ['GCA', 'GCC', 'GCU', 'GCG'],
                'N': ['AAU', 'AAC'],
                'F': ['UUC', 'UUU'],
                'Y': ['UAC', 'UAU'],
                '$': ['UGA', 'UAA', 'UAG'],
                'L': ['UUG', 'CUA', 'CUU', 'CUC', 'CUG', 'UUA'],
                'I': ['AUC', 'AUU', 'AUA'],
                'P': ['CCU', 'CCG', 'CCC', 'CCA'],
                'Q': ['CAG', 'CAA'],
                'M': ['AUG'],
                'S': ['AGC', 'AGU', 'UCG', 'UCU', 'UCC', 'UCA'],
                'E': ['GAA', 'GAG'],
                'C': ['UGC', 'UGU'],
                'K': ['AAA', 'AAG'],
                'R': ['CGA', 'CGG', 'AGG', 'CGC', 'AGA', 'CGU'],
                'W': ['UGG'],
                'T': ['ACU', 'ACC', 'ACG', 'ACA'],
                'V': ['GUA', 'GUU', 'GUC', 'GUG'],
                'D': ['GAU', 'GAC'],
                'G': ['GGU', 'GGA', 'GGG', 'GGC'],
                'H': ['CAU', 'CAC']}

DnaComplimentary = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T': 'A'}

DnaToRna = {'A' : 'A', 'C' : 'C', 'G' : 'G', 'T': 'U'}

RnaToDna = {'A' : 'A', 'C' : 'C', 'G' : 'G', 'U': 'T'}

def GetReversComplimentaryDna(dna):
    result = ""
    for i in range(len(dna)):
        result += DnaComplimentary[dna[i]]

    return result[::-1]

####################        Функции        ####################

# Для выода диагностических сообщейний при дебаге
def printdbg(*args):
    if DEBUG_:
        print(*args)

#  Временный класс, мб вообще удалю или перепишу
class SimpleReport:
    flag = 0
    description = ""
    def __init__(self, fl, descr):
        self.flag = fl
        self.description = descr
    def __repr__(self):
        return [self.flag, self.description]

# Для нашей работы нам наужно сам геном и файл выравнивания с белками
def TotalCheck(_genome, _protein_alignment):

    # Набор функций для проверки по отдельным критериям.
    # Тк мы работаем над небольшими геномами, то кажется можем позволить себе стратегию модульной проверки
    # те отдельной проверки результата по каждому критерию независимо. Каждая функция проверки будет возвращать
    # подмножество из _protein_alignment удовлетворяющее ее критерию.
    # Итоговый результат должен быть объединением всех промежуточных проверок - это и будет список надежных генов.

    def CheckOpenReadingFrame(genome, protein_alignment):
        return SimpleReport(0, "FAILURE")

    def CheckExoneIntroneStructure(genome, protein_alignment):
        return SimpleReport(0, "FAILURE")

    def CheckRibosomalBindingSite(genome, protein_alignment):
        return SimpleReport(0, "FAILURE")

    def CheckCodonUsageBias(genome, protein_alignment):
        return SimpleReport(0, "FAILURE")

    check_stack = [CheckOpenReadingFrame,
                  CheckExoneIntroneStructure,
                  CheckRibosomalBindingSite,
                  CheckCodonUsageBias]

    # В разработке
    report = []
    for function in check_stack:
        report.append(function(_genome, _protein_alignment))

    return report

# {имя_белка: [(start_pos, stop_pos), ...]}
def ReadProteinAlignment(filename = "", file_from = "exonerate"):
    if filename == "":
        Exception("Set the file name for check")

    file_from = file_from.lower()

    file = open(filename, 'r')

    str_file = file.read()

    result = {}
    # В разработке, поэтому пишем парсер спрева только для exonerate'a
    if file_from == "exonerate":
        # Парсим exonerate
        first_token = "Query:"
        second_token = "Target range:"
        offset = 1
        pos = offset
        protein_name = ""
        left_adr = ""
        right_adr = ""
        while True:
            # Получаем имя белка
            pos = str_file.find(first_token, pos)
            if pos == -1:
                break
            pos += len(first_token)
            pos += offset
            protein_name = ""
            while str_file[pos] != '\n':
                protein_name += str_file[pos]
                pos += 1

            # Получаем адреса выравнивания
            left_adr = ""
            right_adr = ""
            pos = str_file.find(second_token, pos)
            if pos == -1:
                break
            pos += len(second_token)
            pos += offset
            # Левый
            while str_file[pos] != ' ':
                left_adr += str_file[pos]
                pos += 1
            pos += len("-> ") + 1
            # Правый
            while str_file[pos] != ' ':
                right_adr += str_file[pos]
                pos += 1

            if protein_name in result.keys():
                result[protein_name].append([int(left_adr), int(right_adr)])
            else:
                result[protein_name] = [int(left_adr), int(right_adr)]
    else:
        Exception("Now available only for the exonerate")

    return result





####################    Основная логика    ####################
def __main__():
    # Парсим аргументы командной строки и вообще ведем себя как порядочная утилита
    # parser = argparse.ArgumentParser()
    #
    # parser.add_argument('-sf', required = True,  help = 'Path to source fasta file')
    # parser.add_argument('-a', required = True, help = 'Path to alignment file')
    # parser.add_argument('-pn', required = True, choices = ["exonerate"] , help = 'Program name what make alignment')
    # parser.add_argument('-d', required = False, help = 'Destination folder for output')
    #
    # parse_result = vars(parser.parse_args())
    #
    # printdbg("Your input:", parse_result)
    #
    # source_fasta   = parse_result["sf"]
    # alignment_file = parse_result["a"]
    # program_name   = parse_result["pn"]
    # destination    = parse_result["d"] if parse_result["d"] else "output"
    #
    # printdbg("Real destination for output:", destination)
    #
    # # Общий формат данных которые мы должны получать из любого выравнивания
    # # имя белка, позиции начала и конца белка в геноме, позиции интронов
    # common_protein_alignment_data = ReadProteinAlignment(alignment_file, program_name)
    # result = TotalCheck(source_fasta, common_protein_alignment_data)
    #
    # try:
    # #   Сохраняем отчет о проверке
    #     output_file = open(destination, "w")
    #
    #     printdbg(result)
    #
    #     output_file.close()
    # except:
    #     print("Can't save a result to " + destination)
    alignment_file = "exonerateout"
    program_name = "exonerate"

    common_protein_alignment_data = ReadProteinAlignment(alignment_file, program_name)
    print(common_protein_alignment_data)

# Запуск
__main__()
