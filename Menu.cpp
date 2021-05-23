#include <iostream>
#include <limits>
#include "Menu.h"

using namespace std;

Menu::Menu(const string &title) {
    this->title = title;
}

unsigned int Menu::pushOption(const string &option) {
    options.push_back(option);
    return options.size() - 1;
}

void Menu::printOptions() {
    for (int i = 0; i < title.size(); i++) {
        cout << "-";
    }
    cout << endl << title << endl;
    for (int i = 0; i < title.size(); i++) {
        cout << "-";
    }
    for (int i = 0; i < options.size(); i++) {
        cout << " " << i + 1 << ": " << options.at(i) << endl;
    }
}

unsigned int Menu::chooseOption() {
    printOptions();
    cout << "Escolha: ";
    unsigned choice;
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cin >> choice;
    while (cin.fail() || cin.peek() != '\n' || choice < 1 || choice > options.size()) {
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Escolha inválida! Insira um nº entre 1 e " << options.size() << ".\n";
        cout << "Escolha: ";
        cin >> choice;
    }
    return choice - 1;
}

const string &Menu::getTitle() const {
    return title;
}

void Menu::setTitle(const string &tit) {
    title = tit;
}

const vector<string> &Menu::getOptions() const {
    return options;
}

void Menu::setOptions(const vector<string> &opts) {
    options = opts;
}
