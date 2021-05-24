#include <iostream>
#include <limits>
#include "Menu.h"
#include "Utils.h"

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
    cout << endl;
    for (int i = 0; i < options.size(); i++) {
        cout << " " << i + 1 << ": " << options.at(i) << endl;
    }
}

unsigned int Menu::chooseOption() {
    printOptions();
    unsigned choice = getUnsigned("Escolha");
    while (choice < 1 || choice > options.size()) {
        cout << "Escolha invalida! Insira um nr entre 1 e " << options.size() << ".\n";
        choice = getUnsigned("Escolha");
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
