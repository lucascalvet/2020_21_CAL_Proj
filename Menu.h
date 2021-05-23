#ifndef INC_2020_21_CAL_PROJ_MENU_H
#define INC_2020_21_CAL_PROJ_MENU_H

#include <string>
#include <vector>

class Menu {
private:
    std::string title;
    std::vector<std::string> options;
public:
    explicit Menu(const std::string &title);

    unsigned pushOption(const std::string& option);

    void printOptions();

    unsigned chooseOption();

    const std::string &getTitle() const;

    void setTitle(const std::string &title);

    const std::vector<std::string> &getOptions() const;

    void setOptions(const std::vector<std::string> &options);
};


#endif //INC_2020_21_CAL_PROJ_MENU_H
