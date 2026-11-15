#!/bin/bash

# 1. Arrancar el agente SSH
eval "$(ssh-agent -s)"

# 2. Cargar tu clave privada (usa el nombre real de tu clave)
ssh-add ~/.ssh/id_ed25519
# o si tu clave se llama id_rsa:
# ssh-add ~/.ssh/id_rsa

# 3. Verificar
ssh-add -l
ssh -T git@github.com