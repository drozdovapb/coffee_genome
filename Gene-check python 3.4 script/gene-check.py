import argparse
# Размер генома кофе в памяти чуть менее 500 мегабайт, можно грузить его туда полностью

# Глобальная переменная, выставляйте в True при дебаге
DEBUG_ = True

# Для выода диагностических сообщейний только при дебаге
def printdbg(*kwargs):
    if DEBUG_:
        print(kwargs)

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


def ReadProteinAlignment(filename = "", file_from = "exonerate"):
    if filename == "":
        Exception("Set the file name for check")

    file_from = file_from.lower()
    # В разработке, поэтому пишем парсер спрева только для exonerate'a

    if file_from == "exonerate":
        pass

# Основная логика скрипта
def __main__():
    # Парсим аргументы командной строки и вообще ведем себя как порядочная утилита
    parser = argparse.ArgumentParser()

    parser.add_argument('-sf', required = True,  help = 'Path to source fasta file')
    parser.add_argument('-a', required = True, help = 'Path to alignment file')
    parser.add_argument('-pn', required = True, choices = ["exonerate"] , help = 'Program name what make alignment')
    parser.add_argument('-d', required = False, help = 'Destination folder for output')

    parse_result = vars(parser.parse_args())

    printdbg("Your input:", parse_result)

    source_fasta   = parse_result["sf"]
    alignment_file = parse_result["a"]
    program_name   = parse_result["pn"]
    destination    = parse_result["d"] if parse_result["d"] else "output"

    printdbg("Real destination for output:", destination)

    # Общий формат данных которые мы должны получать из любого выравнивания
    # имя белка, позиции начала и конца белка в геноме, позиции интронов
    common_protein_alignment_data = ReadProteinAlignment(alignment_file, program_name)
    result = TotalCheck(source_fasta, common_protein_alignment_data)

    try:
    #   Сохраняем отчет о проверке
        output_file = open(destination, "w")

        printdbg(result)

        output_file.close()
    except:
        print("Can't save a result to " + destination)

# Запуск
__main__()