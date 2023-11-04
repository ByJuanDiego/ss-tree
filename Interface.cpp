#include <SFML/Graphics.hpp>
#include <vector>
#include <sstream>
#include <SFML/Network.hpp>

#include "tinyfiledialogs.h"

#include "SStree.h"
#include "CortexAPI.h"

class Button {
public:
    Button(float x, float y, float width, float height, const std::string &text)
            : isHovered(false) {

        shape.setPosition(sf::Vector2f(x, y));
        shape.setSize(sf::Vector2f(width, height));
        shape.setFillColor(sf::Color::White);
        shape.setOutlineThickness(2.0f);
        shape.setOutlineColor(sf::Color::Black);

        font = sf::Font();
        if (!font.loadFromFile("content/Ubuntu-M.ttf")) {
            std::cerr << "No se pudo cargar la fuente!" << std::endl;
        }
        buttonText.setFont(font);
        buttonText.setString(text);
        buttonText.setCharacterSize(14);
        buttonText.setFillColor(sf::Color::Black);
        buttonText.setPosition(
                x + (width - buttonText.getLocalBounds().width) / 2,
                y + (height - buttonText.getLocalBounds().height) / 2
        );
    }

    bool isClicked(sf::Event e) const {
        if (e.type == sf::Event::MouseButtonPressed) {
            if (e.mouseButton.button == sf::Mouse::Left) {
                if (shape.getGlobalBounds().contains(e.mouseButton.x, e.mouseButton.y)) {
                    return true;
                }
            }
        }
        return false;
    }

    void draw(sf::RenderWindow &window) {
        if (shape.getGlobalBounds().contains(sf::Mouse::getPosition(window).x, sf::Mouse::getPosition(window).y)) {
            if (!isHovered) {
                isHovered = true;
                shape.setFillColor(sf::Color(200, 200, 200));
            }
        } else {
            if (isHovered) {
                isHovered = false;
                shape.setFillColor(sf::Color::White);
            }
        }

        window.draw(shape);
        window.draw(buttonText);
    }

    void setSize(float width, float height) {
        shape.setSize(sf::Vector2f(width, height));
        centerText();
    }

    void setPosition(float x, float y) {
        shape.setPosition(sf::Vector2f(x, y));
        centerText();
    }

private:
    sf::RectangleShape shape;
    sf::Text buttonText;
    sf::Font font;
    bool isHovered;

    void centerText() {
        sf::FloatRect textBounds = buttonText.getGlobalBounds();
        float xText = shape.getPosition().x + (shape.getSize().x - textBounds.width) / 2;
        float yText = shape.getPosition().y + (shape.getSize().y - textBounds.height) / 2;
        buttonText.setPosition(xText, yText);
    }
};


class ImageSearchApp {
public:
    ImageSearchApp();
    ~ImageSearchApp() {
        for (sf::Sprite* sprite: resultSprites) {
            delete sprite;
        }
        for (sf::Texture* texture: resultTextures) {
            delete texture;
        }
    }
    void run();

private:
    void init();
    void processEvents();
    void render();

    void loadImage();
    void searchImages();
    void resizeSpriteTo(sf::Sprite &sprite, float width, float height);

private:
    sf::RenderWindow window;
    sf::Texture selectedTexture;
    sf::Sprite selectedSprite;
    std::vector<sf::Texture*> resultTextures;
    std::vector<sf::Sprite*> resultSprites;
    SsTree sstree;
    CortexAPI cortex;
    Button selectButton;
    Button searchButton;
    bool imageSelected = false;
    char const *filepath_of_selected_image = NULL;
};

void ImageSearchApp::init() {
    selectButton.setSize(150, 40);
    selectButton.setPosition(10, 40);
    searchButton.setSize(150, 40);
    searchButton.setPosition(170, 40);
}

ImageSearchApp::ImageSearchApp()
        : window(sf::VideoMode(1200, 800), "Buscador de Imágenes"),
          selectButton(10, 10, 100, 50, "Seleccionar"),
          searchButton(120, 10, 100, 50, "Buscar"),
          sstree(448){
    sstree.loadFromFile("embedding.dat");
    init();
}

void ImageSearchApp::run() {
    while (window.isOpen()) {
        processEvents();
        render();
    }
}

void ImageSearchApp::processEvents() {
    sf::Event event;
    while (window.pollEvent(event)) {
        if (event.type == sf::Event::Closed) {
            window.close();
        }
        if (selectButton.isClicked(event)) {
            loadImage();
        }
        if (searchButton.isClicked(event) && imageSelected) {
            searchImages();
        }
    }
}

void ImageSearchApp::render() {
    window.clear(sf::Color(200, 200, 200));

    resizeSpriteTo(selectedSprite, 400, 400);
    selectedSprite.setPosition(10, 90);

    int xPos  = 420;
    int yPos  = 90;
    int index = 0;

    for (sf::Sprite *sprite : resultSprites) {
        resizeSpriteTo(*sprite, 200, 200);
        sprite->setPosition(xPos, yPos);
        window.draw(*sprite);
        index++;
        xPos += 210;
        if (index % 3 == 0) {
            xPos = 420;
            yPos += 210;
        }
    }

    if (imageSelected) {
        window.draw(selectedSprite);
    }
    selectButton.draw(window);
    searchButton.draw(window);
    window.display();
}

void ImageSearchApp::loadImage() {
    const char* filters[] = {"*.jpg", "*.png", "*.bmp", "*.jpeg"};
    const char* filepath = tinyfd_openFileDialog("Selecciona una imagen", "", 4, filters, NULL, 0);

    if (filepath && selectedTexture.loadFromFile(filepath)) {
        selectedSprite.setTexture(selectedTexture);
        imageSelected = true;
        filepath_of_selected_image = filepath;
    }
}

void ImageSearchApp::searchImages() {
    if (imageSelected) {
        std::vector<NType> imageVec = cortex.postImage(filepath_of_selected_image);
        auto paths = sstree.kNNQuery(Point(imageVec), 6);

        resultTextures.clear();
        resultSprites.clear();

        for (const auto &path : paths) {
            sf::Texture* texture = new sf::Texture();
            if (texture->loadFromFile(path)) {
                resultTextures.push_back(texture);
                sf::Sprite* sprite = new sf::Sprite();
                sprite->setTexture(*texture);
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
