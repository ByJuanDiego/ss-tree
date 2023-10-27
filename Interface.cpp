#include <SFML/Graphics.hpp>
#include <vector>
#include <sstream>
#include <SFML/Network.hpp>

#include "tinyfiledialogs.h"

#include "SStree.h"
#include "CortexAPI.h"

class ImageSearchApp {
public:
    ImageSearchApp();

    void run();

private:
    void processEvents();
    void update();
    void render();

    void loadImage();
    void searchImages();
    void resizeSpriteTo(sf::Sprite &sprite, float width, float height);

private:
    sf::RenderWindow window;
    sf::Texture selectedTexture;
    sf::Sprite selectedSprite;
    std::vector<sf::Texture> resultTextures;
    std::vector<sf::Sprite> resultSprites;
    SsTree sstree;
    CortexAPI cortex;
};

ImageSearchApp::ImageSearchApp() : window(sf::VideoMode(1200, 800), "Buscador de Imágenes"), sstree(448) {
    sstree.loadFromFile("sstree.dat");
}

void ImageSearchApp::run() {
    while (window.isOpen()) {
        processEvents();
        update();
        render();
    }
}

void ImageSearchApp::processEvents() {
    sf::Event event;
    while (window.pollEvent(event)) {
        if (event.type == sf::Event::Closed)
            window.close();
        if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::L) {
            loadImage();
        }
        if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::S) {
            searchImages();
        }
    }
}

void ImageSearchApp::update() {
    // Actualización lógica de tu aplicación
}

void ImageSearchApp::render() {
    window.clear();

    // Ajustar tamaño de la imagen seleccionada y mostrarla
    resizeSpriteTo(selectedSprite, 500, 500);
    window.draw(selectedSprite);

    // Mostrar los 6 resultados en 2 filas y 3 columnas
    int xPos = 550;
    int yPos = 0;
    int index = 0;

    for (const auto &sprite : resultSprites) {
        resizeSpriteTo(resultSprites[index], 250, 250);
        resultSprites[index].setPosition(xPos, yPos);
        window.draw(resultSprites[index]);

        index++;
        xPos += 260;

        if (index % 3 == 0) {
            xPos = 550;
            yPos += 260;
        }
    }

    window.display();
}

void ImageSearchApp::loadImage() {
    const char *filters[] = {"*.jpg", "*.png", "*.bmp"};
    char const *filepath = tinyfd_openFileDialog("Selecciona una imagen", "", 3, filters, NULL, 0);

    if (filepath && selectedTexture.loadFromFile(filepath)) {
        selectedSprite.setTexture(selectedTexture);
        selectedSprite.setPosition(0, 0);
    }
}

void ImageSearchApp::searchImages() {
    const char *filters[] = {"*.jpg", "*.png", "*.bmp"};
    const char *filepath = tinyfd_openFileDialog("Selecciona una imagen", "", 3, filters, NULL, 0);

    if (filepath) {
        std::vector<NType> imageVec = cortex.postImage(filepath);
        auto paths = sstree.kNNQuery(Point(imageVec), 6);

        // Asumo que paths retorna una lista de rutas de imágenes
        resultTextures.clear();
        resultSprites.clear();
        
        for (const auto &path : paths) {
            sf::Texture texture;
            if (texture.loadFromFile(path)) {
                resultTextures.push_back(texture);
                sf::Sprite sprite(texture);
                resultSprites.push_back(sprite);
            }
        }
    }
}

void ImageSearchApp::resizeSpriteTo(sf::Sprite &sprite, float width, float height) {
    float scaleX = width / sprite.getLocalBounds().width;
    float scaleY = height / sprite.getLocalBounds().height;
    sprite.setScale(scaleX, scaleY);
}

int main() {
    ImageSearchApp app;
    app.run();
    return 0;
}
