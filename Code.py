import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

# Initialisation de toutes les données :

# Création du dictionnaire qui contient toutes les informations des voitures

voitures = {
    1: (1760, 5.1, 5.28, 1.95, 1.35, 0.38, 0.3, 0.1),
    2: (1615, 5, 4.51, 1.81, 1.27, 0.29, 0.3, 0.1),
    3: (1498, 5.3, 4.72, 1.88, 1.30, 0.35, 0.3, 0.1),
    4: (1385, 5.2, 4.3, 1.75, 1.23, 0.28, 0.3, 0.1),
    5: (1540, 5.8, 4.6, 1.79, 1.36, 0.34, 0.3, 0.1),
    6: (1600, 5, 4.51, 1.81, 1.48, 0.28, 0.3, 0.1),
}


# Fonction pour accéder + facilement aux caractéristiques

def Get_Attribute(car, attribute):
    match attribute:
        case "Masse":
            return voitures[car][0]
        case "Accélération":
            return voitures[car][1]
        case "Longueur":
            return voitures[car][2]
        case "Largeur":
            return voitures[car][3]
        case "Hauteur":
            return voitures[car][4]
        case "Cx":
            return voitures[car][5]
        case "Cz":
            return voitures[car][6]
        case "µ":
            return voitures[car][7]


def course(car):
    kr = 0.002
    m = Get_Attribute(car, "Masse")
    g = 9.81
    v0 = 0
    angle = 40
    Cx = Get_Attribute(car, "Cx")
    Surface = 3 * 10 ** -4
    p_air = 1.225
    k = 0.5 * Cx * Surface * p_air
    h_looping = 0.23
    r_looping = h_looping / 2
    hauteur_saut = -0.1
    longueur_saut = 0.7
    T = 0.001
    P = 0.001
    temps_total = 0

    # region Pente :
    print("Calcul pour la pente :")
    h_pente = 2

    taille_pente = h_pente / math.sin(angle * math.pi / 180)

    a_sans_frottements = g * math.sin(angle * math.pi / 180)
    a_avec_frottements = g * math.sin(angle * math.pi / 180) - k * g * math.cos(angle * math.pi / 180)

    t_general = np.linspace(0, taille_pente, int(taille_pente * 100))

    vp_sans_frottement = [v0]
    for i in range(1, int(taille_pente * 100)):
        v = math.sqrt(v0 ** 2 + 2 * (i / 100) * a_sans_frottements)
        vp_sans_frottement.append(v)

    vp_avec_frottement = [v0]
    for i in range(1, int(taille_pente * 100)):
        v = math.sqrt(v0 ** 2 + 2 * (i / 100) * a_avec_frottements)
        vp_avec_frottement.append(v)

    plt.plot(t_general, vp_sans_frottement, label="Sans frottements")
    plt.plot(t_general, vp_avec_frottement, label='Avec frottements')
    plt.title("Vitesse de la voiture en fonction de la distance :")
    plt.xlabel("Distance avec le point de départ en m")
    plt.ylabel("Vitesse en m/s")
    plt.legend()
    plt.grid()
    plt.show(block=False)
    plt.pause(5)
    plt.close()

    print("Vitesse à la fin de la pente sans frottements :", vp_sans_frottement[-1])
    print("Vitesse à la fin de la pente avec frottements :", vp_avec_frottement[-1])

    temps_total += taille_pente / vp_avec_frottement[-1]
    # endregion
    # region Looping

    print("Calcul pour le looping :")

    v0_sans_frottements_looping = vp_sans_frottement[-1]
    v0_avec_frottements_looping = vp_avec_frottement[-1]

    v0_min_pour_resussir = math.sqrt(5 * g * r_looping)
    print("Vitesse minimale pour réussir le looping :", v0_min_pour_resussir)
    hauteur_min_pour_reussir = (v0_min_pour_resussir ** 2 * math.sin(angle * math.pi / 180)) / (2 * a_sans_frottements)
    print("Hauteur minimale pour réussir le looping :", hauteur_min_pour_reussir)

    if v0_sans_frottements_looping > v0_min_pour_resussir and v0_avec_frottements_looping > v0_min_pour_resussir:

        nb_elements = 360

        vl_sans_frottements = [v0_sans_frottements_looping]
        teta = np.linspace(0, 360, nb_elements)

        for i in range(0, len(teta) - 1):
            vl_sans_frottements.append(r_looping * math.sqrt((v0_sans_frottements_looping / r_looping) ** 2 - (
                    2 * g * (1 - math.cos(teta[i] * math.pi / 180))) / r_looping))

        plt.plot(teta, vl_sans_frottements, label="Sans frottements")
        plt.title("Vitesse de la voiture en fonction de l'ange téta")
        plt.xlabel("Evolution de l'angle téta")
        plt.ylabel("Vitesse en m/s")
        plt.legend()
        plt.grid()
        plt.show(block=False)
        plt.pause(5)
        plt.close()

        v_sortie_looping_sf = vl_sans_frottements[-1]

        # Looping avec frottements :

        x0 = 0
        S_0 = (x0, v0_avec_frottements_looping)

        def dSdx(x, S):
            x, v = S
            return [v, (-m * g * math.sin(x) + (r_looping * v ** 2 * m - m * g * math.cos(x)) * kr - k * (
                    (r_looping * v) ** 2)) / (m * r_looping)]

        teta = np.linspace(0, 0.689, nb_elements)
        sol = odeint(dSdx, y0=S_0, t=teta, tfirst=True)

        x_sol = sol.T[0]
        v_sol = sol.T[1]

        plt.plot(teta, v_sol, color="blue", lw=2, label="Avec frottements")
        plt.title("Vitesse de la voiture dans le looping en fonction du temps :")
        plt.xlabel("Temps (s)")
        plt.ylabel("Vitesse en m/s")
        plt.legend()
        plt.grid()
        plt.show(block=False)
        plt.pause(5)
        plt.close()

        v_sortie_looping_f = v_sol[-1]
        print("Vitesse en sortie de looping sans frottements :", v_sortie_looping_sf)
        print("Vitesse en sortie de looping avec frottements :", v_sortie_looping_f)

        distance = 2 * np.pi * 6
        temps_total += distance / v_sortie_looping_f

        # endregion
        # region Saut:
        # Partie saut sans frottements :
        print("Calcul pour le saut :")

        # Déterminons t :
        t = math.sqrt((hauteur_saut * -2) / g)

        # Déterminons la vitesse minimale pour réussir le saut :
        v_min = longueur_saut / t
        print("Vitesse minimale pour réussir le saut :", v_min)

        v0 = v_sortie_looping_sf
        duree = np.linspace(0, t)

        # Vitesse de la voiture :
        v_mobile = []
        for i in range(0, len(duree)):
            v_mobile.append(math.sqrt((-g * duree[i]) ** 2 + v0 ** 2))

        # Position de la voiture :
        x = [0]
        for i in range(len(duree)):
            x.append(v0 * duree[i])

        y = [0]
        for i in range(len(duree)):
            y.append(-1 / 2 * g * duree[i] ** 2)

        # Affichage des graphiques :

        plt.plot(duree, v_mobile)
        plt.title("Vitesse de la voiture sans frottements :")
        plt.xlabel("Distance avec le point de départ en m")
        plt.ylabel("Vitesse en m/s")
        plt.grid()
        plt.show(block=False)
        plt.pause(5)
        plt.close()

        f = plt.figure()
        ax = f.add_subplot(111)
        plt.axhline(y=-0.1)
        plt.axvline(x=0.7)
        plt.plot(x, y, color="orange")
        plt.title("Trajectoire de la voiture sans frottements :")
        plt.xlabel("Distance point de départ en m")
        plt.ylabel("Hauteur en m")
        plt.grid()
        plt.show(block=False)
        plt.pause(5)
        plt.close()

        # Partie saut avec frottements :

        v0x = v_sortie_looping_f  # vitesse initiale en x
        v0y = 0  # vitesse initiale en y
        t0 = 0
        nb_echantillons = 100

        def vpp(vp, t):
            vx = vp[2]
            vy = vp[3]

            ax = -p_air / (2 * m) * np.sqrt(vx ** 2 + vy ** 2) * (T * vx + P * vy)
            ay = -p_air / (2 * m) * np.sqrt(vx ** 2 + vy ** 2) * (T * vy - P * vx) - g

            return [vx, vy, ax, ay]

        PositionVitesse0 = [0, 0, v0x, v0y]
        t = np.linspace(t0, 0.3, nb_echantillons)

        PositionVitesse = odeint(vpp, PositionVitesse0, t)

        f = plt.figure()
        ax = f.add_subplot(111)
        plt.axhline(y=-0.1)
        plt.axvline(x=0.7)
        plt.plot(PositionVitesse[:, 0], PositionVitesse[:, 1], color="orange")
        plt.title("Trajectoire de la voiture avec frottements :")
        plt.xlabel("Distance du point de départ (m)")
        plt.ylabel("Hauteur (m)")
        plt.grid()
        plt.show(block=False)
        plt.pause(5)
        plt.close()
        temps_total += t[-1]
        # endregion:
        # region Piste

        # Paramètres de la voiture
        masse = Get_Attribute(car, "Masse")  # kg
        acceleration_moyenne = Get_Attribute(car, "Accélération")  # m/s^2
        largeur = Get_Attribute(car, "Largeur")  # m
        hauteur = Get_Attribute(car, "Hauteur")  # m
        cx = Get_Attribute(car, "Cx")
        mu = Get_Attribute(car, "µ")
        g = 9.81  # m/s^2, accélération due à la gravité

        # Conditions initiales
        vitesse_initiale = 0.0  # m/s
        position_initiale = 0.0  # m

        # Fonction représentant le système d'équations différentielles
        def equation_mouvement(y, t):
            vitesse, position = y
            am = acceleration_moyenne
            lambda_val = 0.5 * cx * 1.225 * largeur * hauteur / masse
            dvdt = am - mu * g - lambda_val * vitesse ** 2 / (2 * masse)
            dydt = [dvdt, vitesse]
            return dydt

        # Temps de simulation
        temps = np.linspace(0, 10, 1000)

        # Résolution des équations différentielles
        resultats = odeint(equation_mouvement, [vitesse_initiale, position_initiale], temps)

        # Extraction de la vitesse et de la position
        vitesse = resultats[:, 0]
        position = resultats[:, 1]

        # Recherche du temps de traversée
        temps_de_traversee = temps[np.argmax(position >= 10)]

        # Affichage des résultats
        plt.plot(temps, position, label='Position')
        plt.xlim(0, 3)
        plt.ylim(0, 11)
        plt.xlabel('Temps (s)')
        plt.ylabel('Position (m)')
        plt.title('Simulation du mouvement de la voiture')
        plt.legend()
        plt.show()
        plt.pause(5)
        plt.close()

        temps_total += temps_de_traversee
        print("Calcul du temps total en cours : ")
        print("le temps total de la course est : " + str(temps_total))
        # endregion
    else:
        print("Le passage du looping n'est pas possible !")


course(int(input("Entrez le numéro du modèle de la voiture : ")))