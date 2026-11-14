## Usuarios Linux

En Linux (Ubuntu) podemos encontrar la consola en el icono llamado **Terminal**, que ya viene instalada en la computadora. Puedes ubicarla al  hacer click en  el icono **Show Apps**

<img src="FIGURES/01_ShowApps.png" alt="Show Apss" style="zoom:50%;" />



En el icono de búsqueda introduce el texto **Terminal** y haz click en el icono para desplegarla

![Type Terminal](FIGURES/02_LookForTerminal.png)

### Atajo

También puedes abrir la terminal directamente con el comando **Ctrl + Alt + T**



### Conexion al servidor

Una vez abierta la **Terminal**, deberás conectarte al clúster ingresando tu usuario y contraseña mediante `ssh`:

```
ssh usuario@cluster.lancis.ecologia.unam.mx
```

![Comando de conexión SSH](FIGURES/03_ssh_connection.png)



Ingresa la contraseña asociada a tu usuario. Para ello puedes copiarla  (Ctrl + C) y pegarla (Ctrl + shift + V)

**PRECACIÓN** Si ingresas erróneamente la contrasena tres veces consecutivas, tu usuario sera bloqueado.

![Ingreso de contraseña en Terminal](FIGURES/04_enter_password.png)

Al ingresar al clúster recibirás el mensaje de bienvenida junto con un breve manual que podrás explorar usando el comando `more`. Con esto estarás list@ para trabajar durante el taller:

![Pantalla de bienvenida del clúster](FIGURES/05_access_cluster.png)
