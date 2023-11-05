# ss-tree

Dependencies

```ssh
sudo pacman -S git-lfs
sudo pacman -S nlohmann-json
sudo pacman -S sfml
sudo pacman -S curl
```

Instructions to set the workspace

```ssh 
git clone https://github.com/ByJuanDiego/ss-tree.git

cd ss-tree
git lfs install
git lfs pull
unzip img.zip
```

Then, compile the interface

```ssh
g++ -std=c++17 Interface.cpp CortexAPI.cpp SStree.cpp tinyfiledialogs.c -o Interface -lsfml-graphics -lsfml-window -lsfml-system -lcurl
```

Finally, call the executable

```ssh
./Interface
```

Search the K-NNs!

![image](https://github.com/ByJuanDiego/ss-tree/assets/79115974/323a2c25-2300-4ece-80eb-65fe263224d0)
