# Videre arbeid og ting som bør sjekkes

Status 26.10.2026. Punktene i del 1 er ting i dagens kode eller dokumentasjon som gir feil eller usikre svar. Del 2 er forslag til verifisering mot litteratur. Del 3 er forslag til videreutvikling. Tall merket «sjekket» er kjørt med dagens kode.

## Anbefalt rekkefølge

1. **Rett `main.py`** (1.1 og 1.3). Den gir feil svar i dag, og retting er rask.
2. **Lag en pytest-suite** (3.1), så resten kan endres trygt.
3. **Gjør steg 1 konsistent med steg 2** (1.4 og 3.2). Det påvirker alle tallene fra rørsimuleringene.
4. **Sammenlign med publiserte løsninger** (del 2).
5. **Utvid med lettere rør, sykluser og egenvekt** (3.3–3.6).

---

## 1. Bør sjekkes eller rettes

### 1.1 `main.py` gir feil bæreevne for stripe som standard
- **Hva:** Standardverdien for ruhet er r = 1,0, og `PrandtlSolver` er bare riktig for r = 0. Med r = 1 gir den q = cᵤ·(1 + π/2) = 2,57·cᵤ. Riktig svar for ru stripe er også cᵤ·(π + 2) = 5,14·cᵤ (sjekket).
- **Forslag:** Sett standardverdien til r = 0, og gi feilmelding eller advarsel for r > 0 til 1.2 er løst.

### 1.2 `PrandtlSolver` for ru stripe (r > 0)
- **Hva:** Region I har konstant spenning med ψ_A = π/2 − ½·asin r. Det er bare konsistent for glatt fundament (står i docstringen i `yield_function.py`). `MohrCoulombYield` gir `NotImplementedError` for r > 0.
- **Forslag:** Lag en egen mekanisme for ru fundament (Hill-mekanisme, der fundamentflaten er en glidelinje eller en stiv kile ligger under), eller avvis r > 0.

### 1.3 Ru skive: vertikal last med skjærretningen fra sideveis skyv
- **Hva:** `DiskPrandtlSolver` bruker som standard skjærretningen fra en horisontal forskyvning. Da er ikke `collapse_load_Q()` konsistent med spenningsfeltet. Ved p = 0,25 og r = 0,5 (sjekket):

  | | Q [kN/m] |
  |---|---|
  | `collapse_load_Q()`, standardretning | 3,32 |
  | 2·F_v fra kontaktspenningene, standardretning | 2,82 |
  | `vertical_shear=True` | 4,23 |

- **Brukes slik i dag:**
  - `main.py` (skive med r > 0) viser Q fra standardretningen.
  - `plot_disk_Fv_vs_p.py`, `plot_disk_Fv_Fh_vs_p.py` og `plot_disk_forces_vs_z.py` bruker standardretningen for r = 0,5 og 1, men figurene kalles «vertical resistance».
- **Forslag:**
  - Bruk `vertical_shear=True` der lasten er vertikal.
  - Gå gjennom figurene i `plots/` som er laget slik.
  - Gjør `collapse_load_Q()` konsistent med retningen, eller gi metoden et navn som sier hva den er.
  - Sjekk også antagelsen F_v(p) = F_v(1 − p) for «extraction» i de samme skriptene. Det er ikke opplagt at opptrekk er symmetrisk med nedtrykk.

### 1.4 Steg 1 overvurderer bæreevnen til det ru røret
- **Hva:** Spenningsløsningen i steg 1 (`DiskPrandtlSolver(vertical_shear=True)`) gir V = 2,80 kN/m ved likevekt. Den ensidige mekanismen fra steg 2 (en øvre grense) bærer bare 2,69–2,70 kN/m ved H = 0, og det holder seg med finere nett. Den sanne bæreevnen er derfor høyst 2,70, og steg 1 ligger minst ~4 % for høyt.
- **Følger:**
  - Røret synker for lite i steg 1.
  - Steg 2 starter med H = −0,29 kN/m i stedet for 0.
  - Retningen endres brått de første millimeterne.
- **Forslag:** Se 3.2.

### 1.5 Hastighetsfeltet for vertikal skive er ikke kinematisk tillatt
- **Hva:** d < 0 i 30–60 % av det plastiske området (p = 0,1–0,5, også for r = 0; sjekket). Spenningsløsningen er altså ikke en komplett løsning, og det er ikke bevist at den er en nedre grense.
- **Størrelse:** Den øvre grensen fra det samme feltet ligger ~10 % over (4,01 mot 3,64 ved p = 0,25).
- **Forslag:** Finn en kinematisk tillatt mekanisme for vertikal belastning. Den gir en ekte øvre grense, og sammen med spenningsløsningen et intervall for den sanne verdien. Dette henger sammen med 1.4.

### 1.6 Nγ i eksempel 4 i manualen: faktor 2 og feil sammenligning
- **Hva:** I koden er B halv bredde, og `q_over_cB` = Q_half/(c·B) = midlere trykk/c. Manualen trekker ut Nγ fra Q/(cB) = Nc + ½(γB/c)·Nγ med halv bredde. Det gir det dobbelte av den vanlige Nγ, som er definert med hel bredde.
- **Følger:** Verdiene 19–25 tilsvarer derfor ~10–13. Sammenligningen med Vesic (22,4) er ikke gyldig. Vesic er dessuten en empirisk formel; de eksakte verdiene er lavere. Fra hukommelsen (sjekk kilden) er Nγ = 7,65 for glatt og 14,75 for ru fundament ved φ = 30° (Martin 2005).
- **Forslag:** Rett definisjonen i eksemplet, og sammenlign med eksakte verdier. Superposisjonen Nc + Nγ er heller ikke eksakt. Jeg oppdaterte bare de forventede tallene i dette eksemplet da manualen ble skrevet om, ikke definisjonen.

### 1.7 β hopper ved friksjon som snur
- **Hva:** β hopper ca. 1° når vendepunktet φ_s flytter seg mellom nodene. Med halvert steglengde endres w og H med ~2 % (146,8 mot 144,3 mm og 0,892 mot 0,874 kN/m).
- **Forslag:** Legg noder nøyaktig på overgangssonens ender φ_s ± ω_r, slik `OneSidedDiskSolver` allerede gjør for θ₀ i sonen med parallelle linjer. Da følger nettet φ_s kontinuerlig.

### 1.8 Overflaten foran røret ved forskyvning tilbake
- **Hva:** Ved snuen tilpasses overflaten foran (grøftas vegg og bermen bak) med én rett linje, α = −24°. Grøfta er ujevn, så det er en grov forenkling akkurat der H er minst sikker.
- **Forslag:** Plott overflaten mot den tilpassede linja de første 10 mm av steg 3. Vurder en knekt overflate eller kortere tilpasningslengde.

### 1.9 Konvergens som ikke er dokumentert
- **`PipePenetration`:** følsomhet for `dw`, `dx` og `n_disk`. Volumfeilen er sjekket (~1e−17 m²), men ikke hvordan synkningen avhenger av steglengden.
- **`PipeLateral`:** hele forløpet med finere nett (`n_disk`, `n_fan`). Bare V ved H = 0 er sjekket med `n_disk = 90`.
- **Forslag:** Lag en liten konvergenstabell for hver og legg den i manualen.

### 1.10 Allerede sjekket og i orden
- **Mohr-Coulomb-skiva** nærmer seg Tresca-skiva når φ → 0: 0,3–0,4 % avvik ved φ = 0,1° og ~4 % ved 1° (sjekket for p = 0,1–0,45 og r = 0 og 0,5).
- **Virtuelt arbeid** stemmer med spenningsløsningen innen 0,02–0,5 % (sjekket).
- **Alle eksemplene i manualen** gir de tallene som står der (sjekket 26.10.2026).

---

## 2. Verifisering mot litteratur

Referansene er skrevet fra hukommelsen. Sjekk tittel, år og tall i kildene før de brukes.

- **Martin & White (2012)**, *Limit analysis of the undrained bearing capacity of offshore pipelines*, Géotechnique 62(9). Øvre og nedre grenser for V, H og V–H-kurver for glatt og ru rør ved ulike w/D. Bruk dem til å kontrollere V(d) i steg 1 (punkt 1.4–1.5) og H ved start av steg 2.
- **Randolph & White (2008)**, *Upper-bound yield envelopes for pipelines at shallow embedment in clay*, Géotechnique 58(4). Mekanismer med kile og skje for sideveis last. Sammenlign med den ensidige mekanismen i `PipeLateral`.
- **Aubeny, Shi & Murff (2005)**, *Collapse loads for a cylinder embedded in trench in cohesive soil*, Int. J. Geomech. 5(4). Vertikal bæreevne mot innbygging.
- **Merifield, White & Randolph (2009)**, *Effect of surface heave on response of partially embedded pipelines on clay*, J. Geotech. Geoenviron. Eng. 135(6). Effekten av heving på vertikal bæreevne. Sammenlign med `PipePenetration`, der hevingen reduserer synkningen med 37–43 %.
- **Wang, White & Randolph (2010)**, *Large-deformation finite element analysis of pipe penetration and large-amplitude lateral displacement*, Can. Geotech. J. 47(8). Forløp ved stor sideforskyvning. Tunge rør (V/V_max nær 1) dykker og lette rør stiger; vårt rør dykker, som forventet for V/V_max = 1. Sammenlign forløpet og H(u) kvalitativt.
- **Martin (2005)**, *Exact bearing capacity calculations using the method of characteristics* (IACMAG Torino). Eksakte Nc og Nγ, til punkt 1.6.

---

## 3. Forslag til videreutvikling

### 3.1 Pytest-suite (lite arbeid, stor nytte)
Gjør dagens kontroller om til tester med `assert`:
- glatt Prandtl: 2 + π;
- `MohrCoulombPrandtlSolver` mot analytisk Nc;
- `DiskPrandtlSolver` mot lukket form for Q;
- Mohr-Coulomb-skiva → Tresca-skiva når φ → 0;
- `energy_balance()`: virtuelt arbeid ≈ spenningsløsning;
- volumbevaring i `PipePenetration`;
- speilsymmetri i `PipeLateral.reversed`;
- eksemplene i manualen.

Hold de trege rørsimuleringene utenfor, eller merk dem som `slow`.

### 3.2 Steg 1 som er konsistent med steg 2
Regn synkningen med en øvre grense i stedet for spenningsløsningen. Mekanismen kan være en symmetrisk variant av den med friksjon som snur og en kile som henger fast under rørbunnen, eller den beste av ensidig og tosidig mekanisme. Da starter steg 2 i likevekt (H = 0). Sammenlign V(d) med Martin & White (2012).

### 3.3 V–H-kurve ved fast geometri
Regn H_max for flere V ved samme innbygging, og plott kurven. Det er en direkte sammenligning med litteraturen og et nyttig designverktøy.

### 3.4 Lette rør (V/V_max < 1)
I praksis synker et rør under installasjonslast og skyves så sideveis med lavere vekt. Det kan legges inn som:
- `PipePenetration` til V_install;
- `PipeLateral` med W < V_install.

Lette rør skal stige og eventuelt bryte ut (breakout). Det er den viktigste situasjonen for lateral knekking (lateral buckling).

### 3.5 Flere sykluser
`reversed()` kan i prinsippet brukes flere ganger, men det er bare testet for én snu. Kjør flere fram-og-tilbake-sykluser og sjekk hvordan bermene bygges opp og hvordan H endrer seg per syklus.

### 3.6 Egenvekt og skjærstyrke som øker med dybden
Det viktigste som mangler for reell offshore leire:
- **Egenvekt γ':** oppdrift og bermenes vekt. Det enkleste er et oppdriftsledd (som hos Merifield mfl.); det grundigste er et Tresca-nett med tyngde.
- **Skjærstyrke sᵤ = sᵤ₀ + k·z:** krever et nytt Riemann-steg, siden Hencky-relasjonene ikke lenger gjelder.

### 3.7 Heving med skrå overflate i steg 1
`PipePenetration` regner med flat overflate på kontaktnivå. `OneSidedDiskSolver` kan allerede håndtere en skrå overflate (α); bruk det samme i steg 1 så stegene blir konsistente.

### 3.8 Gjenskapbart dashbord
Dataene til dashbordet ble eksportert med skript og pickle-filer i en midlertidig mappe for økten, og de kan forsvinne. Flytt eksportskriptet inn i repoet (f.eks. `export_dashboard.py`), så dashbordet kan bygges på nytt fra `run_lateral.py`.

### 3.9 Ytelse
`run_lateral.py` tar ~3 min (friksjon som snur) og ~7 min (fri ω). Mulige forbedringer:
- Kjør β-søket parallelt (`multiprocessing`), siden Mac-en har flere kjerner.
- Gjenbruk geometrien til nettet mellom nabosteg.

### 3.10 Rydding
- **ruff:** De trygge rettingene fra `ruff check --fix` er gjort (3.11.2026, samme utskrift fra alle skriptene før og etter). Det er 39 funn igjen: 35 × C408 (`dict(...)` i stedet for `{...}`), 2 × BLE001 (`except Exception`) og 2 × RUF046. De krever `--unsafe-fixes` eller retting for hånd; ta det i én egen commit, så ikke ekte endringer drukner. `ruff format` er kjørt på alle filene (3.11.2026).
- **Gamle skript:** vurder å fjerne eller slå sammen `dbg.py`, `test_disk2.py` og liknende.
- **`solver.py`:** ~2000 linjer; kan deles opp (stripe, skive, Mohr-Coulomb).
- **Lokale endringer:** bestem hva som skal skje med de endrede `plots/mc_gravity_*.png` og `manual.pdf`.

### 3.11 Manualen
- Legg til en figur for sideforskyvningen (nett, hastighet og forløp). I dag har den bare tabeller.
- Ta inn resultatene fra 1.4–1.6 når de er avklart.
