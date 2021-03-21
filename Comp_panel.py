import math


class Material:
    def __init__(self, sigma_ult, sigma_yield, E, density):
        self.sigma_ult = sigma_ult
        self.sigma_yield = sigma_yield
        self.E = E
        self.density = density


class Stringer:
    def __init__(self, a, t, material, length):
        self.a = a
        self.t = t
        self.material = material
        self.length = length
        self.area = a * t + (a - t) * t
        self.I = (1 / 12) * a * pow(t, 3) + pow((0.5 * t), 2) * t * a + (1 / 12) * t * pow((a - t), 3) + (
                a - t) * t * pow((0.5 * a + t), 2)
        self.mass = self.area * length * material.density
        self.hole_mass = t * pow(0.0032 / 2, 2) * math.pi * material.density


class Skin:
    def __init__(self, a, b, t, material):
        self.a = a
        self.b = b
        self.t = t
        self.material = material
        self.area = t * b
        self.I = (1 / 12) * pow(t, 3) * b
        self.mass = a * b * t * material.density
        self.hole_mass = t * pow(0.0032 / 2, 2) * math.pi * material.density


class Rivet:
    def __init__(self, diametre, lenght, material, F_max_sheer):
        self.d = diametre
        self.lenght = lenght
        self.material = material
        self.F_max_sheer = F_max_sheer
        self.mass = lenght * pow((diametre / 2), 2) * material.density * math.pi


class Panel:
    def __init__(self, config, rivet_per_stringer, profiles, Skin, req_buckle_force, req_ult_force):
        self.config = config
        self.rivet_per_stringer = rivet_per_stringer
        self.config = config
        self.profiles = profiles
        self.Skin = Skin
        self.req_buckle_force = req_buckle_force
        self.req_ult_force = req_ult_force
        self.rivet = riv_short
        if self.Skin.t > 0.00119 and self.profiles[1] != 0:  # should be skin.t + stringer
            self.rivet = riv_long
        self.rivet_count = sum(config) * self.rivet_per_stringer
        self.total_mass = self.rivet_count * (self.rivet.mass - self.Skin.hole_mass) + self.Skin.mass
        for x in range(len(self.profiles)):
            self.total_mass += (config[x] * profiles[x].mass - self.rivet_per_stringer * profiles[x].hole_mass)

    def ultimate_check(self):
        total_area = self.Skin.area
        for x in range(len(self.profiles)):
            total_area += self.config[x] * self.profiles[x].area
        return total_area * self.Skin.material.sigma_yield > self.req_ult_force  # rivets are irrelevant

    def buckle_check(self):
        if sum(config) < 5:
            return False
        rivet_spacing = (self.profiles[0].length - 0.01) / (
                self.rivet_per_stringer + 1)  # check if +1 is correct
        y_bar_individual = []
        I_total_individual = []
        boundry_condition = 4
        pop_rivet_c = 2.1
        buckle_a_difference = 0.0  # What's that? It shows up here and no other value is ever assigned to it
        stringer_spacing = self.Skin.b / (sum(config) - 1)
        a_by_b = (self.Skin.a - buckle_a_difference) / stringer_spacing
        K_c = 6.3
        if a_by_b < 6:
            K_c = 0.0114 * pow(a_by_b, 4) - 0.239 * pow(a_by_b, 3) + 1.84 * pow(a_by_b, 2) - 6.1744 * a_by_b + 14.017
        total_area = self.Skin.area
        stringer_area = 0
        for x in range(len(self.profiles)):
            stringer_area += self.config[x] * self.profiles[x].area
        total_area += stringer_area
        # checkpoint
        for x in range(len(self.profiles)):
            y_bar_individual.append(y_bar_cal(self.Skin, self.profiles[x]))
        y_bar_skin = y_bar_individual[0][0]  # same skin for all iterations
        Ay_bar_stringers = 0
        for x in range(len(self.config)):
            Ay_bar_stringers += (self.config[x] * y_bar_individual[x][1] * self.profiles[x].area)
        y_bar = (y_bar_skin * self.Skin.area + Ay_bar_stringers) / total_area
        # checkpoint II
        for x in range(len(self.profiles)):
            I_total_individual.append(
                I_total_cal(self.Skin, self.profiles[x], y_bar, y_bar_cal(self.Skin, self.profiles[x])))
        I_total = I_total_individual[0][0]
        for x in range(len(self.config)):
            I_total += self.config[x] * I_total_individual[x][1]
        buckle_force_column = buckle_force_column_cal(self.Skin.a - buckle_a_difference, self.Skin.material.E,
                                                      boundry_condition, I_total)
        applied_buckle_stress = (self.req_buckle_force) / (total_area)
        sigma_crit_sheet = K_c * self.Skin.material.E * pow((self.Skin.t / stringer_spacing), 2)
        sigma_ir = 0.9 * pop_rivet_c * self.Skin.material.E * pow((self.Skin.t / rivet_spacing), 2)
        sigma_buckle = buckle_force_column / total_area
        global sheet
        global inner
        global column
        if sigma_crit_sheet > applied_buckle_stress:
            sheet += 1
        elif sigma_ir > applied_buckle_stress:
            inner += 1
        elif sigma_buckle > applied_buckle_stress:
            column += 1

        return sigma_crit_sheet > applied_buckle_stress and sigma_ir > applied_buckle_stress and sigma_buckle > applied_buckle_stress


def y_bar_cal(Skin, Stringer):
    SkinT = Skin.t
    SkinB = Skin.b  # It doesn't change skin y_bar, why do we want it here?
    StrT = Stringer.t
    StrA = Stringer.a
    # I have no idea what's going on here
    # y_bar_skin = (((SkinT * SkinB) / (SkinT * SkinB + StrT * (StrA + StrA - StrT))) * (0.5 * SkinT))
    # y_bar_stringer = (((StrT * (StrA + StrA - StrT)) / (SkinT * SkinB + StrT * (StrA + StrA - StrT))) * (
    #        (((StrA - StrT) * StrT) / (StrT * (StrA + StrA - StrT))) * (0.5 * StrA + 0.5 * StrT) + (
    #        (StrT * StrA) / (StrT * (StrA + StrA - StrT))) * (0.5 * StrT)))
    y_bar_skin = SkinT * 0.5
    y_bar_stringer = SkinT + (StrT * 0.5 * (StrA - StrT) * StrA + 0.5 * StrA * StrA * StrT) / (
            StrA * StrT + (StrA - StrT) * StrT)
    return [y_bar_skin, y_bar_stringer]


def I_total_cal(Skin, Stringer, y_bar, y_bar_arr):
    SkinT = Skin.t
    SkinB = Skin.b
    StrT = Stringer.t
    StrA = Stringer.a
    ybsk = y_bar_arr[0]
    ybstr = y_bar_arr[1] - SkinT
    yb = y_bar
    # I_total_skin = (((1 / 12) * SkinB * (pow(SkinT, 3))) + ((pow((yb - 0.5 * SkinT), 2)) * SkinB * SkinT))
    # I_total_stringer = (((1 / 12) * StrT * (pow((StrA - StrT), 3))) + (
    # (pow(((StrA * 0.5 + 0.5 * StrT) - yb), 2)) * SkinT * (StrA - SkinT)) + (
    #  ((1 / 12) * StrA * (pow(StrT, 3)) + (pow((yb - 0.5 * StrT), 2)) * StrT * StrA)))
    I_total_skin = (1 / 12) * SkinB * (pow(SkinT, 3)) + Skin.area * pow((yb - SkinT * 0.5), 2)
    I_stringer = (1 / 12) * (StrA - StrT) * pow(StrT, 3) + (1 / 12) * StrT * pow(StrA, 3) + (StrA - StrT) * StrT * pow(
        ybstr - StrT * 0.5, 2) + StrT * StrA * pow(ybstr - StrA * 0.5, 2)  # Stringer about its y_bar
    I_total_stringer = I_stringer + Stringer.area * pow(yb - ybstr, 2)  # Stringer about panel y_bar
    return [I_total_skin, I_total_stringer]


def buckle_force_column_cal(a, E, boundry_condition, I_total):
    return (boundry_condition * E * I_total * pow(math.pi, 2)) / (pow(a, 2))


def mass_sort(Panel):
    return Panel.total_mass


inner = 0
sheet = 0
column = 0
aluminium = Material(127000000, 100000000, 63500000000, 2780)
steel = Material(400000000, 25000000, 200000000000, 7850)
L1 = Stringer(0.02, 0.0015, aluminium, 0.495)  # changed lenght of the stringrs
L2 = Stringer(0.02, 0.002, aluminium, 0.495)
L3 = Stringer(0.015, 0.001, aluminium, 0.495)
L4 = Stringer(0.015, 0.0015, aluminium, 0.495)
F1 = Skin(0.495, 0.4, 0.0008, aluminium)
F2 = Skin(0.495, 0.4, 0.0012, aluminium)
riv_short = Rivet(0.0032, 0.006, steel, 1060)
riv_long = Rivet(0.0032, 0.0104, steel, 1060)
config = [7, 0, 0, 0]
profiles = [L1, L2, L3, L4]
Fbuckling = 34500*1.2
Fult = 60000

print(y_bar_cal(F1, L4))

# panel = Panel(config, 12, profiles, F1, 34500, 60000)

possible_config = []
mass_guard = 1

for x in range(0, 12, 1):
    for y in range(0, 12, 1):
        for z in range(0, 12, 1):
            for w in range(0, 12, 1):
                for s in range(2, 20, 1):
                    config = [x, y, z, w]
                    panel = Panel(config, s, profiles, F2, Fbuckling, Fult)
                    if panel.total_mass < mass_guard:
                        if panel.buckle_check() and panel.ultimate_check():
                            possible_config.append(panel)
                            mass_guard = panel.total_mass
for x in range(0, 12, 1):
    for y in range(0, 12, 1):
        for z in range(0, 12, 1):
            for w in range(0, 12, 1):
                for s in range(2, 20, 1):
                    config = [x, y, z, w]
                    panel = Panel(config, s, profiles, F1, Fbuckling, Fult)
                    if panel.total_mass < mass_guard:
                        if panel.buckle_check() and panel.ultimate_check():
                            possible_config.append(panel)
                            mass_guard = panel.total_mass

possible_config.sort(key=mass_sort)

for x in range(len(possible_config)):
    print(str(possible_config[x].config) + ' ' + str(possible_config[x].rivet_per_stringer) + ' ' + str(
        possible_config[x].total_mass) + ' ' + str(possible_config[x].Skin.t))
print([sheet, inner, column])