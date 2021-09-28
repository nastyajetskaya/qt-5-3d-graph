#include "method.h"
#include "window.h"
#include <QApplication>
#include <QMessageBox>

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    Window scene;

    if (scene.parse_command_line(argc, argv)) {
        QMessageBox::warning(0, "Wrong input arguments!",
                             "Wrong input arguments!");
        return -1;
    }

    scene.setWindowTitle("Approx3D");
    scene.resize(500, 500);
    scene.showMaximized();
    return app.exec();
}
