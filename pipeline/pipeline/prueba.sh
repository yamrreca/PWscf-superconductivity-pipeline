#!/bin/bash

read
if [[ "$REPLY" == "y" || "$REPLY" == "n" ]]; then
	echo "Hola mundo"
else
	echo "Respuesta inválida"
fi
