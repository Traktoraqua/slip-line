# Plan: dynamisk modell for røret

Status 13.11.2026. Planen innfører én ting av gangen. Hvert steg er ferdig først når dokumentasjonen er oppdatert og kontrollene er kjørt. Ingen steg startes før forrige steg er godkjent.

## Status

| Steg | Innhold | Status |
|---|---|---|
| 0 | Plan, notasjon og valg | ferdig |
| 1 | Flyteflaten (H–V) for én tilstand | ferdig |
| 2 | Bevegelseslikningen på fast flyteflate | ferdig |
| 3 | Geometrien følger bevegelsen (forover) | ferdig, venter på godkjenning |
| 4 | Fram og tilbake | ikke startet |
| 5 | Harmonisk bølgelast og løft | ikke startet |
| 6 | Jordas treghet | ikke startet |
| 7 | Hastighetsavhengig cᵤ | ikke startet |

---

## Arbeidsmåte for hvert steg

Dokumentasjonen kommer først og er en del av steget, ikke noe som gjøres til slutt.

1. **Teori først.** Skriv likningene og antakelsene i manualen før koden. Lukkede formler med utledning skrives i `pipe_formulas.tex`.
2. **Kode** i `slip_line/`, med docstrings og engelske navn.
3. **Kontroller** med tall vi vet svaret på. Der det passer, som pytest-tester.
4. **Figur** fra et `plot_*.py`-skript i repoet, lagret i `plots/`.
5. **Oppdater dokumentene:**
   - `manual.tex`: teori, API og et eksempel med forventet utskrift. Kjør eksemplet på nytt.
   - `README.md` og `CLAUDE.md`.
   - `VIDERE_ARBEID.md`, hvis et punkt der løses eller endres.
   - Statustabellen i denne fila.
6. **Bygg og sjekk:**
   - Manualen bygges med 0 feil og ingen nye overfulle linjer, og figurene sjekkes visuelt.
   - `ruff format --check` og `ruff check` er rene.
7. **Commit** i små steg. Ikke push.
8. **Oppsummering** til Trygve, og klarsignal før neste steg.

## Hvor dokumentasjonen står

| Hva | Hvor |
|---|---|
| Teori, API, eksempler og figurer | `manual.tex`, ny seksjon «Dynamic response of the pipe» etter «Lateral push of a pipe» (engelsk, som resten av manualen) |
| Lukkede formler med utledning og kontroll | `pipe_formulas.tex` (norsk) |
| Status, valg og rekkefølge | denne fila |
| Kjente svakheter og usikre resultater | `VIDERE_ARBEID.md` |
| Kort oversikt og kommandoer | `README.md` og `CLAUDE.md` |

## Notasjon

- **u, w:** horisontal forskyvning (mot +x) og vertikal forskyvning (nedover) av rørsenteret.
- **δ:** bevegelsesretningen under horisontalen, med n = (cos δ, sin δ). Akse y peker nedover.
- **D(δ):** dissipasjonen per fartsenhet når røret beveger seg i retning δ, fra mekanismen.
- **(H, V):** kraften fra røret på jorda. Flyteflaten er alle (H, V) jorda tåler: H·cos δ + V·sin δ ≤ D(δ) for alle δ.
- **Fx(t), F_L(t):** ytre last på røret, horisontalt og løft (oppover).
- **W:** rørets neddykkede vekt.
- **m:** rørets masse med innhold og belegg.
- **m_a,x, m_a,z:** tilleggsmasse fra vannet, horisontalt og vertikalt.
- **M:** massen i bevegelseslikningen, M = diag(m + m_a,x, m + m_a,z). Fra steg 6 kommer jordas masse i tillegg.
- **J(δ):** ∫|v/v₀|² dA over den plastiske sonen. Jordas masse er ρ_sat·J.
- **T, T_n:** bølgeperioden og egenperioden til røret på jorda, T_n = 2π·√(M/k), der k = dH/du langs skyvet ved V = W. T_n avhenger mye av cᵤ, se steg 3.

## Valg som er tatt

- **Frihetsgrader og materialmodell.** To frihetsgrader (u, w). Røret er et stivt legeme, og jorda er stiv-plastisk: ingen elastisk fase og ingen demping utover det plastiske arbeidet.
- **Jordmotstand.** Den kommer fra flyteflaten til den sideveis mekanismen (øvre grense), med normalitet: bevegelsen står vinkelrett på flyteflaten. Der flyteflaten har et hjørne, kan røret bevege seg i alle retninger mellom normalene på hver side.
- **Q fra synkningen brukes ikke.** Når røret har flyttet seg sideveis, åpnes kontakten bak røret. Ved V = W gir Q dessuten samme svar som flyteflaten uten Q.
  - Svakhet: flyteflaten er for sterk når røret beveger seg bratt nedover. Ved start tåler den V = 3,27 kN/m ved H = 0 med cᵤ = 5 kPa og 3,78 kN/m med cᵤ = 1 kPa, mot Q = W = 2,80 kN/m. Den er altså 17 % og 35 % for sterk. Dette har bare betydning når røret lander etter å ha lettet.
- **Friksjonsmodell.** Lukket form med jevn friksjon (`PipeLateralClosedForm`) kommer først, fordi den er rask: én flyteflate med 722 retninger tar 11–14 ms (målt). Friksjon som snur kan komme senere.
- **Tidsintegrasjon.** Implisitt steg med projeksjon på flyteflaten, som «closest point projection» i plastisitetsteori.
  - Prøvekraften er F* = F_ytre + M·v/Δt.
  - Ligger F* innenfor flyteflaten, står røret stille.
  - Ellers er jordkraften punktet på flyteflaten nærmest F*, målt med M⁻¹. Den nye farten blir v = Δt·M⁻¹·(F* − F_jord).
  - Samme regel tar hånd om ro, start, stopp, hjørner og at røret letter.

---

## Steg 1 – Flyteflaten (H–V) for én tilstand

Dette dekker også punkt 3.3 i `VIDERE_ARBEID.md`.

- **Mål:** jordmotstanden i alle retninger for en gitt geometri.
- **Lages:** `slip_line/envelope.py` med `YieldEnvelope`. Den inneholder:
  - D(δ) forover, fra tilstandens mekanisme.
  - D(δ) bakover, fra speilbildet (`reversed`).
  - D = 0 der røret slipper jorda.
  - Skjæringen av halvplanene, som gir flyteflaten som en mangekant.
  - `support(δ)`, som gir kraften for en bevegelsesretning, og `contains(F)`.
  - At den også kan lages direkte fra en mangekant, til tester.
  - Koden flyttes fra skriptene i scratchpad (`envelope.py`, `plot_envelope.py`).
- **Kontroll:**
  - Ved V = W gir flyteflaten samme H og δ som `mechanism()`, for flere tilstander: start, 50 mm, 200 mm og etter skyv tilbake. Ved start er det 0,155 kN/m og 17°.
  - Ved start er forover og bakover like (0,155 og 0,156 kN/m).
  - Flyteflaten går gjennom origo (ingen strekk) og er konveks.
- **Dokumentasjon:**
  - Manualen: underseksjon «Yield envelope», med øvre grense per retning, halvplan, normalitet, hjørner, og hvorfor Q ikke er med. API og eksempel.
  - Figur `plots/yield_envelope.png` fra `plot_yield_envelope.py`: flyteflaten ved start og etter 50 og 200 mm, så veksten og skjevheten synes.
  - `VIDERE_ARBEID.md`: merk 3.3 som løst.
- **Resultat (10.11.2026):**
  - `slip_line/envelope.py` med `YieldEnvelope`. I tillegg til det som er nevnt over har den `push_force(V)`, skyvet forover med styrt forskyvning, slik `mechanism()` regner, og `H_range`, `V_range`, `directions` og `active`.
  - Kontrollene ligger i `tests/test_yield_envelope.py` (8 tester, 1 s). De er skrevet for pytest, men kjøres som skript til pytest er lagt inn.
  - Manualen: ny seksjon 8 «Dynamic response of the pipe» med teori, resultater, API og eksempel. `README.md`, `CLAUDE.md` og `VIDERE_ARBEID.md` (3.3 løst, 3.1 og ny 1.10) er oppdatert.
  - **Avvik fra planen:**
    - Skyvet styrer forskyvningen, så røret får ikke gå bakover. Derfor er det `push_force(W)` som skal være lik `mechanism()`, ikke den største H på flyteflaten. De to er like så lenge halvplanene bakover ikke skjærer flyteflaten ved V = W.
    - Med cᵤ = 5 kPa har røret løftet seg og mistet kontakten bak seg allerede etter ~8 mm, og fra 30 mm når flyteflaten ikke V = W. Figuren viser derfor cᵤ = 5 kPa ved 0, 4 og 10 mm, og kjøringen med cᵤ = 1 kPa ved 0, 50 og 200 mm. Kontrollen ved 50 mm, 200 mm og etter skyv tilbake bruker også cᵤ = 1 kPa.
  - **Funn:**
    - Med cᵤ = 5 kPa følger den kvasistatiske kjøringen et lokalt minimum fra 1,25 mm. H er der opptil 43 % for høy (`VIDERE_ARBEID.md` 1.10).
    - Når røret mangler kontakt bak seg, blir flyteflaten ensidig: H ≥ 0, og ved H = 0 bærer den ingenting. Glipa bak røret er likevel bare brøkdeler av en millimeter.

## Steg 2 – Bevegelseslikningen på fast flyteflate

- **Mål:** røret som stivt legeme med masse på stiv-plastisk jord. Geometrien holdes fast, så bare bevegelseslikningen er ny.
- **Lages:**
  - `slip_line/dynamics.py` med projeksjonen og tidsstegene.
  - Lasten som funksjoner av tid, Fx(t) og F_L(t).
  - Tidsserier med t, u, w, u̇, ẇ, H, V og δ.
- **Kontroll (pytest).** Det krever `uv add --dev pytest`, som installeres lokalt i prosjektet og ikke globalt. `tests/test_yield_envelope.py` fra steg 1 er allerede skrevet for pytest. pytest må bare lete i `tests/`, siden `test_*.py` i rotmappa er skript.
  - Last innenfor flyteflaten: røret står stille.
  - Glideblokk (Newmark): flyteflaten er et rektangel |H| ≤ H₀, 0 ≤ V ≤ V₀ med V₀ > W, og Fx er harmonisk med amplitude F₀ > H₀. Ved V = W er kanten loddrett, så røret glir vannrett, og projeksjonen gir (m + m_a,x)·ü = Fx − H₀. Tidspunktet røret stopper finnes av én likning (rotfinning), og forskyvningen per halvperiode er da en formel. Numerikken skal treffe innen 0,1 %.
    - Ikke H ≤ μ·V: med normalitet beveger røret seg da vinkelrett på kanten, altså oppover i vinkelen atan μ. Kanten går gjennom origo, så D = 0. Det er ingen glideblokk, og ingenting dissiperes.
  - Løft større enn W: røret letter uten jordkraft.
  - Energibalanse: arbeidet fra lastene er lik dissipasjonen pluss den kinetiske energien, innen feilen fra Δt.
  - Konvergens med Δt.
- **Dokumentasjon:**
  - Manualen: underseksjon «Equations of motion», med likningene og projeksjonen forklart med en figur. API og eksempel.
  - Glideblokkformelen som kontroll.
- **Data:** M. Som eksempel brukes 1900 kg/m (1090 kg/m fra W, pluss C_a = 1).
- **Resultat (13.11.2026):**
  - `slip_line/dynamics.py` med `PipeDynamics` (`step`, `run`, `energy`), `closest_point`, `DynamicState` og `pipe_mass`. Ett tidssteg tar ~10 µs på flyteflaten til røret (32 hjørner).
  - Kontrollene ligger i `tests/test_dynamics.py` (12 tester, 0,2 s). pytest er lagt inn: `uv run pytest` kjørte da alle 20 testene på ~1 s.
  - Manualen: underseksjon 8.4 «Equations of motion» med teori, figur og kontroller, og 8.5 med API og eksempel. Figur `plots/projection.png` fra `plot_projection.py`. `README.md`, `CLAUDE.md` og `VIDERE_ARBEID.md` (3.1 delvis gjort) er oppdatert.
  - **Avvik fra planen:**
    - `pipe_mass(W, D, C_a=1)` er lagt til, så C_a kan endres senere. Massen er m = W/g + ρ_w·πD²/4 = 285 + 805 = 1090 kg/m: vekten i vann pluss oppdriften, med D = 1 m og sjøvann.
    - Energibalansen går opp eksakt når den numeriske dempingen i det implisitte steget regnes som et eget ledd. Det leddet halveres når Δt halveres.
    - Figuren viser også glideblokken mot lukket form.
  - **Funn:**
    - Glideblokken treffer lukket form innen 1,5·10⁻⁶ med Δt = T/2000. Formelen per halvperiode gjelder for F₀/H₀ ≤ 1,862; over det snur blokken uten å stoppe.
    - Feilen er av første orden midt i en glidning (0,7 % i u med Δt = T/2000), men mye mindre når blokken har stoppet.
    - På den faste flyteflaten ved start (cᵤ = 5 kPa) flytter en last på 0,2 kN/m (T = 10 s) røret 58 mm fram og tilbake og 35 mm ned per periode. Det er langt utenfor der flyteflaten gjelder: skyvet er 0,398 kN/m allerede etter 4 mm. Tall for virkelige laster krever steg 3.

## Steg 3 – Geometrien følger bevegelsen (forover)

- **Mål:** røret graver seg ned og bygger berm mens det beveger seg, og flyteflaten vokser.
- **Lages:**
  - `advance` for en vilkårlig bevegelsesretning. I dag flytter den bare langs mekanismens β.
  - Ny flyteflate når røret har flyttet seg Δs, ikke for hvert tidssteg.
  - At røret kan lette og lande igjen.
- **Kontroll:**
  - Langsom last (T ≫ T_n) gir samme forløp som den kvasistatiske kraftstyrte kjøringen. Bruk samme lastform som Fx = 0,3·(1 − cos πt), strukket til f.eks. 100 s.
    - Fra steg 1: den gamle kjøringen følger et lokalt minimum, og den lar ikke røret gå bakover. Sammenlign derfor med `push_force` på hele flyteflaten, og bruk cᵤ = 1 kPa, der den gamle kjøringen er riktig til ~190 mm.
    - Bare den stigende halvdelen kan sammenlignes. Røret trenger en holdekraft nesten like stor som skyvet: med cᵤ = 1 kPa er H_min = 0,35 mot H_max = 0,43 kN/m ved 5 mm, og 0,48 mot 0,60 kN/m ved 25 mm. Når lasten faller fra toppen på 0,6 kN/m, går røret bakover og ned mot glipa bak seg (δ = 152° ved 20 mm) allerede ved Fx ≈ 0,47 kN/m. Da er u̇ < 0, og det krever den speilede rammen fra steg 4. Kjøringene i steg 3 stopper derfor når u̇ ville blitt negativ.
  - Røret mister kontakten bak seg når det løfter seg (fra ~8 mm med cᵤ = 5 kPa). Glipa er brøkdeler av en millimeter. Δs må være så liten at røret lander i glipa i stedet for å synke gjennom den, eller flyteflaten må lages på nytt når kontakten bak endres. Med cᵤ = 1 kPa går røret mot glipa bak seg når lasten faller (se over), og den landingen hører til steg 4.
  - Den stigende halvdelen på 1 s (lasten uten strekk) viser hvor mye dynamikken betyr. Med M = 1900 kg/m og k fra `push_force` langs skyvet er T_n ≈ 2,0–3,2 s for cᵤ = 1 kPa (k ≈ 7,5–20 kN/m² de første 25 mm) og 0,7–1,7 s for cᵤ = 5 kPa (k ≈ 27–136 kN/m² de første 4 mm). k er størst den første millimeteren.
  - Volumet bevares som i dag.
  - Energibalanse, og konvergens med Δt og Δs.
- **Dokumentasjon:**
  - Manualen: underseksjon om geometrien i tid, med figur av u(t), w(t) og lastens vei gjennom flyteflaten.
  - Nytt skript `run_dynamic.py`, beskrevet i `README.md` og `CLAUDE.md`.
- **Resultat (13.11.2026):**
  - `PipeLateral.move(state, du, dw, flow)` flytter røret og jorda for en vilkårlig bevegelse, også en blanding av to mekanismer i et hjørne av flyteflaten. `advance` er nå `move` langs tilstandens mekanisme, og `run_lateral_closed_form.py` og `plot_yield_envelope.py` gir nøyaktig samme utskrift som før.
  - `PipeMotion` i `slip_line/dynamics.py` kjører tidsstegene fra steg 2 på flyteflaten til den aktuelle geometrien. Når røret har flyttet seg Δs (0,25 mm), flyttes jorda med mekanismene til de kantene farten har fulgt, og det lages en ny flyteflate. Uten kontakt er jordkraften null, og geometrien oppdateres med en gang røret treffer jorda igjen. Kjøringen stopper før en bevegelse som trenger en speilet mekanisme (`stop_reason = "backward"`).
  - Kontrollene ligger i `tests/test_motion.py` (11 tester, ~10 s). `uv run pytest` kjører alle 32 på ~12 s.
  - Manualen: 8.6 (teori), 8.7 (resultater, figur og kontroller) og 8.8 (API og eksempel). Figur `plots/dynamic_push.png` fra `plot_dynamic.py`, og nytt skript `run_dynamic.py` (~4 s). `README.md`, `CLAUDE.md` og `VIDERE_ARBEID.md` (1.10) er oppdatert.
  - **Avvik fra planen:**
    - Landingen i glipa bak røret med cᵤ = 5 kPa er flyttet til steg 4. Glipa er bare noen få µm, men V = W ligger da utenfor flyteflaten, og å synke tilbake ned i glipa er en blanding av en mekanisme forover og en speilet mekanisme bakover. Kjøringen stopper der, ved u = 1,5 mm.
    - En landing rett ned i det symmetriske hullet etter at røret har lettet helt trenger også mekanismer på begge sider, og stopper derfor kjøringen. Jorda røret når ned i det siste tidssteget før det treffer, kuttes uten flyt, så volumet bevares ikke helt ved en landing.
    - Testen av rask last antok at røret flytter seg kortere enn ved langsom last. Det stemte ikke (se funn), og testen sjekker nå forsinkelsen i stedet.
  - **Funn:**
    - Langsom last (stiger over 50 s) følger det kvasistatiske skyvet innen 1 % i H og 0,05 mm i w, og stopper ved u = 24,45 mm. Røret går i små rykk: etter hver Δs er den nye flyteflaten litt sterkere, og røret venter til lasten tar igjen. Når lasten faller, stopper kjøringen ved holdekraften H_min = 0,479 kN/m (u̇ < 0 trenger steg 4).
    - Rask last (stiger over 1 s, T_n ≈ 2–3 s): røret begynner ved samme last, men massen gjør at det henger etter. Det har størst fart (42 mm/s) etter lasttoppen, bremses av jorda og stopper når lasten har falt til 0,18 kN/m. Det ender på u = 24,97 mm, 2 % lenger enn ved langsom last. V varierer fra 2,76 til 2,88 kN/m.
    - Δs = 0,5, 0,25 og 0,125 mm gir u = 25,19, 24,97 og 24,84 mm (omtrent første orden). Halvert Δt endrer 0,01 mm.
    - Med cᵤ = 5 kPa følger den dynamiske modellen det laveste minimumet: røret løfter seg fra u ≈ 1,3 mm (δ = −2°). Den kvasistatiske kjøringen fant ikke dette (`VIDERE_ARBEID.md` 1.10).

## Steg 4 – Fram og tilbake

Dette dekker også punkt 3.5 i `VIDERE_ARBEID.md`.

- **Mål:** lasten skifter retning, og røret kan gå begge veier mange ganger.
- **Lages:** bytte mellom rammen forover og den speilede rammen bakover når u̇ skifter fortegn.
- **Kontroll:**
  - Symmetrisk last gir speilet svar.
  - Langsom last gir dagens skyv fram og tilbake (steg 3 i `run_lateral_closed_form.py`).
  - Når lasten faller, lander røret i glipa bak seg i stedet for å synke gjennom den. Med cᵤ = 1 kPa skjer det allerede ved Fx ≈ 0,47 kN/m etter 20 mm (se steg 3).
  - Med cᵤ = 5 kPa slipper røret sålen bak seg når det løfter seg (fra u ≈ 1,5 mm), og synker så tilbake ned i glipa med en blanding av en mekanisme forover og en bakover (flyttet fra steg 3). Δs må være så liten at røret lander i glipa i stedet for å synke gjennom den. Sammenlign med `push_force` på flyteflaten.
  - Landing rett ned i hullet etter at røret har lettet helt (begge sider).
  - Bermene bygges opp på begge sider.
- **Dokumentasjon:**
  - Manualen: om vendinger, med figur av geometrien etter noen perioder.
  - `VIDERE_ARBEID.md`: merk 3.5.

## Steg 5 – Harmonisk bølgelast og løft

- **Mål:** bølgelast over hele stormer, og hvor langt røret flytter seg per periode.
- **Lages:**
  - En enkel lastmodell med parametere:
    - Fx(t) = F₀·sin(2πt/T), eventuelt med strøm.
    - F_L(t) ∝ U², som topper to ganger per periode, samtidig med draget.
  - Hopp over tiden røret ligger i ro. Lasten er kjent, så tidspunktet den går ut av flyteflaten kan finnes direkte.
  - Utskrift per periode: forskyvningen, og om røret låser seg eller vandrer.
- **Kontroll:**
  - Hoppingen gir samme svar som fulle tidssteg.
  - T ≫ T_n gir kvasistatisk svar.
  - Skalering: samme F₀/H_maks og T/T_n gir samme u/D.
- **Dokumentasjon:**
  - Manualen: lastmodellen og et eksempel med en storm.
  - Figur av u(t) over mange perioder og lastens vei i flyteflaten.
  - `README.md`.
- **Data fra Trygve:** F₀, T, strøm og løftamplitude.

## Steg 6 – Jordas treghet

- **Mål:** ta med massen til jorda som beveger seg. I eksemplene var den 15–40 % av den totale massen.
- **Lages:**
  - J(δ), først numerisk fra nettet (`InclinedDiskVelocityField`) og så på lukket form for jevn friksjon.
  - J tas med i M, først som ρ_sat·J(δ) for den aktuelle retningen.
- **Kontroll:**
  - Lukket form mot summen over nettets trekanter, med konvergens når n_disk øker, som for arealet.
  - Tallene fra tidligere: J/R² = 0,98 ved P/R = 0,4, r = 0,8, θ = 10° og δ = 30°.
  - Kjør steg 5 på nytt og vis hvor mye svaret endres.
- **Dokumentasjon:**
  - `pipe_formulas.tex`: ny seksjon etter arealet, med utledning.
  - Manualen: om massen.
  - `CLAUDE.md`.
- **Data fra Trygve:** ρ_sat.

## Steg 7 – Hastighetsavhengig cᵤ

- **Mål:** cᵤ øker med tøyningshastigheten, typisk ~10 % per tidobling.
- **Lages:** cᵤ = cᵤ,ref·(1 + μ·log₁₀(γ̇/γ̇_ref)), der γ̇ ≈ v delt på størrelsen av den plastiske sonen.
- **Kontroll:**
  - μ = 0 gir svaret fra steg 6.
  - Vis effekten på forskyvningen per periode.
- **Dokumentasjon:** manualen.
- **Data fra Trygve:** μ og hastigheten cᵤ er målt ved.

---

## Senere (valgfritt)

- **Et ledd ∝ ρv².** Det kommer fra ny jord som krysser den ytre α-linja og får fart brått. Det betyr noe først ved ~1 m/s når cᵤ = 5 kPa.
- **Jordas vekt γ′**, se punkt 3.6 i `VIDERE_ARBEID.md`. Ved cᵤ = 5 kPa er γ′D/cᵤ ≈ 1,1.
- **Vertikal mekanisme for skjev overflate**, for når røret lander etter å ha lettet.
- **Friksjon som snur på lukket form**, med φ_s = π/4 − δ.
- **Vannkrefter fra vannhastigheten** (Morison eller DNV-RP-F109) i stedet for gitte krefter.
- **Stabilitetskurver:** forskyvning mot F₀/H_maks og T/T_n.

## Data og valg som trengs fra Trygve

| Når | Hva | Forslag hvis ikke gitt |
|---|---|---|
| Steg 2 | Rørmasse m og C_a (horisontalt og vertikalt) | 1090 kg/m fra W, og C_a = 1 (avklart 13.11.2026: greit for nå, C_a er en input i `pipe_mass`). C_a = 1 er verdien i fritt vann. Nær bunnen gir DNV-RP-F105 rundt 2,3 for et rør som ligger på bunnen (ikke kontrollert her), og mindre når røret er delvis gravd ned. |
| Steg 2 | Ok å ta inn pytest (`uv add --dev pytest`) | ja (avklart 13.11.2026) |
| Steg 5 | Bølgelast: F₀, T, strøm og løft | T = 10 s, F₀ rundt største motstand |
| Steg 6 | ρ_sat | 1600 kg/m³ |
| Steg 7 | μ og hastigheten cᵤ er målt ved | μ = 0,1 |
