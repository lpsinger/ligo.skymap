.. highlight:: sh

Easter Egg: List Most Probable Constellations (`ligo-skymap-constellations`)
============================================================================

.. argparse::
    :module: ligo.skymap.tool.ligo_skymap_constellations
    :func: parser

To list the most probable constellations for GW170817::

    $ curl -OL https://dcc.ligo.org/public/0146/G1701985/001/bayestar.fits.gz
    $ ligo-skymap-constellations bayestar.fits.gz

.. testcode::
   :hide:

   from ligo.skymap.too.ligo_skymap_constellations import main
   from astropy.utils.data import download_file
   filename = download_file('https://dcc.ligo.org/public/0146/G1701985/001/bayestar.fits.gz', cache=True)
   main([filename])

.. testoutput::

   0.49372269756706605	Virgo
   0.31327241872626527	Corvus
   0.1806012308409065	Hydra
   0.012273063139920255	Centaurus
   7.33302173347542e-05	Lupus
   2.0124397034626024e-05	Taurus
   1.607777623406115e-05	Circinus
   1.1847507101258453e-05	Leo Minor
   4.801076976240868e-06	Apus
   1.3453783636882628e-06	Auriga
   1.1268223149776268e-06	Octans
   1.0061694481116907e-06	Triangulum Australe
   5.674617238901326e-07	Lynx
   3.3689975342661287e-07	Leo
   2.202701360022312e-08	Cetus
   3.4858172839407564e-09	Eridanus
   2.6653058899772474e-10	Tucana
   2.0926237030447548e-10	Ursa Major
   1.523794259929362e-11	Indus
   1.1852612764294879e-11	Coma Berenices
   3.7753942225276276e-12	Perseus
   7.233182922297495e-14	Hydrus
   1.8928276524970434e-15	Aries
   1.806934635454502e-17	Pavo
   1.0647849818860385e-17	Fornax
   1.3310266267130931e-18	Norma
   4.6320972661082874e-23	Cancer
   4.0130685715837047e-23	Gemini
   6.790291479923671e-25	Ara
   3.7764700790357345e-25	Musca
   2.639996865100319e-28	Phoenix
   5.780653404806961e-31	Crater
   1.808312000338703e-31	Chamaleon
   7.113935406971093e-35	Sculptor
   7.794034729999757e-40	Canis Major
   3.002828370601988e-40	Orion
   2.8836795827752243e-40	Vulpecula
   1.2296371584055717e-40	Horologium
   1.123897670645498e-40	Sagitta
   9.996341027424896e-41	Libra
   8.540285308176656e-41	Lyra
   3.641656312609553e-41	Hercules
   3.4677979351969684e-41	Pisces
   3.2810345132622196e-41	Puppis
   3.0667887873006756e-41	Aquila
   2.2835674607055147e-41	Monoceros
   7.913013196981244e-42	Cygnus
   6.929615805137575e-42	Ophiucus
   5.375566642109065e-42	Boötes
   7.55551462044855e-43	Serpens
   5.407878461019846e-44	Delphinus
   1.1364378339902471e-44	Columba
   7.819373074697878e-45	Canes Venatici
   4.2105730126666225e-45	Lepus
   1.9531486764615627e-45	Scutum
   1.4708796040109254e-45	Pyxis
   1.2704687538073003e-45	Camelopardalis
   5.360929060573291e-46	Sagittarius
   2.075911655160898e-46	Canis Minor
   8.481094210479841e-48	Aquarius
   7.705001075439435e-48	Mensa
   5.103078925278023e-48	Capricornus
   1.5173542236819764e-49	Microscopium
   7.560358027520032e-50	Vela
   6.946478828114356e-51	Draco
   1.5475583515525288e-51	Corona Borealis
   3.985453510500382e-52	Equuleus
   3.939003447554475e-52	Crux
   2.0358475122510653e-52	Pictor
   1.2276558370479365e-52	Telescopium
   5.423426832655179e-53	Scorpius
   3.0014643319394047e-53	Triangulum
   1.0712303738745325e-53	Pegasus
   2.5770407860180895e-54	Carina
   1.8380418856152116e-54	Corona Australis
   6.490109434373655e-55	Caelum
   2.582893153639923e-55	Grus
   6.508267657854811e-56	Cepheus
   6.978538846136061e-57	Andromeda
   1.3501808258168752e-57	Sextans
   8.473375548551851e-58	Ursa Minor
   6.233211200719527e-59	Pisces Austrinus
   3.756361457050049e-59	Cassiopeia
   1.3653977982548696e-60	Dorado
   1.1435513819465776e-60	Volans
   6.268628976830198e-64	Reticulum
   1.7967538027310186e-68	Lacerta
   4.6830227320478515e-85	Antlia
