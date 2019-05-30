--
-- PostgreSQL database dump
--

-- Dumped from database version 10.7
-- Dumped by pg_dump version 10.7

-- Started on 2019-05-23 19:10:43

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET client_min_messages = warning;
SET row_security = off;

--
-- TOC entry 2935 (class 0 OID 16565)
-- Dependencies: 203
-- Data for Name: representation; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.representation (id, name) FROM stdin;
0	EMDB contour
1	EMDB contour + 1/3 core
2	EMDB contour + 2/3 core
3	EMDB contour + 1/3 + 2/3 core
4	EMDB contour + 1 std dev
\.


--
-- TOC entry 2944 (class 0 OID 16757)
-- Dependencies: 212
-- Data for Name: user_role; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.user_role (id, role) FROM stdin;
1	User
2	Admin
\.


--
-- TOC entry 2945 (class 0 OID 16808)
-- Dependencies: 213
-- Data for Name: user; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public."user" (id, name, email, role) FROM stdin;
113214863598857275093	Jose Salas	joseantonio122@gmail.com	2
0	Default	1	1
1132148763	Idris Craft	Faint@gmail.com	1
11321644863	Mollie Walmsley	Rustic@gmail.com	1
111322314863	Leena Herman	Miscreant@gmail.com	1
11321448463	Zuzanna Ellison	Dynamic@gmail.com	1
111321245863	Carmel OMoore	user2@gmail.com	1
11323148563	Lennie Reese	Attractive@gmail.com	1
11232514863	Ember Bean	Recondite@gmail.com	1
1134214862	Reanna Dixon	Dangerous@gmail.com	1
1332144864	Elif Skinner	Shiny@gmail.com	1
1413221454865	Leroy Hurley	Dirty@gmail.com	1
19133214866	Desiree Jeffery	Bashful@gmail.com	1
\.


--
-- TOC entry 2937 (class 0 OID 16578)
-- Dependencies: 205
-- Data for Name: volume_filter; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.volume_filter (id, name) FROM stdin;
0	on
1	On
\.


--
-- TOC entry 2947 (class 0 OID 17006)
-- Dependencies: 215
-- Data for Name: benchmark_history; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.benchmark_history (id, date_time, ip, user_id, representation_id, volume_filter_id, top_results, emd_list) FROM stdin;
1	2019-05-24 00:51:21.132	::1	\N	3	1	20	["0002","5502","0017"]
2	2019-05-24 00:51:45.332	::1	\N	0	0	20	["1010","1884","5502"]
\.


--
-- TOC entry 2929 (class 0 OID 16510)
-- Dependencies: 197
-- Data for Name: map_information; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.map_information (id, file_information, data_type, num_columns, num_rows, num_sections, origin_col, origin_row, origin_sec, limit_col, limit_row, limit_sec, spacing_col, spacing_row, spacing_sec, cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma, axis_order_fast, axis_order_medium, axis_order_slow, minimum, maximum, average, std, space_group_number, details, pixel_x, pixel_y, pixel_z, countour_level, annotation_details) FROM stdin;
1	{"format": "CCP4", "sizeKb": 67109, "type": "map", "file": "emd_0001.map.gz"}	Image stored as Reals	256	256	256	0	0	0	255	255	255	256	256	256	{"units": "A", "value": 1.0}	{"units": "A", "value": 1.0}	{"units": "A", "value": 1.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.10788548000000001	0.22688343999999999	1.2941896e-005	0.0077065425	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.017999999999999999	\N
4	{"format": "CCP4", "sizeKb": 67109, "type": "map", "file": "emd_0001.map.gz"}	Image stored as Reals	256	256	256	0	0	0	255	255	255	256	256	256	{"units": "A", "value": 1.0}	{"units": "A", "value": 1.0}	{"units": "A", "value": 1.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.10788548000000001	0.22688343999999999	1.2941896e-005	0.0077065425	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.017999999999999999	\N
5	{"format": "CCP4", "sizeKb": 67109, "type": "map", "file": "emd_0001.map.gz"}	Image stored as Reals	256	256	256	0	0	0	255	255	255	256	256	256	{"units": "A", "value": 1.0}	{"units": "A", "value": 1.0}	{"units": "A", "value": 1.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.10788548000000001	0.22688343999999999	1.2941896e-005	0.0077065425	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.017999999999999999	\N
6	{"format": "CCP4", "sizeKb": 67109, "type": "map", "file": "emd_0002.map.gz"}	Image stored as Reals	256	256	256	0	0	0	255	255	255	256	256	256	{"units": "A", "value": 276.48}	{"units": "A", "value": 276.48}	{"units": "A", "value": 276.48}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.081265143999999997	0.19144611	0.00053819930000000005	0.0061218822000000004	1	\N	{"units": "A", "value": 1.08}	{"units": "A", "value": 1.08}	{"units": "A", "value": 1.08}	0.016	\N
7	{"format": "CCP4", "sizeKb": 202613, "type": "map", "file": "emd_0004.map.gz"}	Image stored as Reals	370	370	370	0	0	0	369	369	369	370	370	370	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-1.2225216999999999	2.4226390000000002	0.0082020170000000007	0.080526760000000003	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.33000000000000002	Modelling map
8	{"format": "CCP4", "sizeKb": 202613, "type": "map", "file": "emd_0005.map.gz"}	Image stored as Reals	370	370	370	0	0	0	369	369	369	370	370	370	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-1.0410119	2.3215709000000002	0.0087016790000000004	0.076787889999999998	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.33000000000000002	Saccharomyces cerevisiae Class1 (III2IV)
9	{"format": "CCP4", "sizeKb": 202613, "type": "map", "file": "emd_0006.map.gz"}	Image stored as Reals	370	370	370	0	0	0	369	369	369	370	370	370	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.57250069999999997	1.5357721	0.0071697747000000001	0.061136416999999998	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.22	Saccharomyces cerevisiae supercomplex class2 (III2IV2)
10	{"format": "CCP4", "sizeKb": 32001, "type": "map", "file": "emd_0007.map.gz"}	Image stored as Reals	200	200	200	0	0	0	199	199	199	200	200	200	{"units": "A", "value": 252.0}	{"units": "A", "value": 252.0}	{"units": "A", "value": 252.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.026832027000000001	0.10347184500000001	0.0019379637	0.0087864729999999995	1	\N	{"units": "A", "value": 1.26}	{"units": "A", "value": 1.26}	{"units": "A", "value": 1.26}	0.042999999999999997	\N
11	{"format": "CCP4", "sizeKb": 364501, "type": "map", "file": "emd_0008.map.gz"}	Image stored as Reals	450	450	450	0	0	0	449	449	449	450	450	450	{"units": "A", "value": 495.0}	{"units": "A", "value": 495.0}	{"units": "A", "value": 495.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	Z	Y	X	-11.851832	30.160060000000001	4.4800000000000003e-012	1	1	\N	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	6.5	\N
12	{"format": "CCP4", "sizeKb": 32001, "type": "map", "file": "emd_0009.map.gz"}	Image stored as Reals	200	200	200	0	0	0	199	199	199	200	200	200	{"units": "A", "value": 220.0}	{"units": "A", "value": 220.0}	{"units": "A", "value": 220.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	Z	Y	X	-9.0799559999999992	16.808064999999999	-2.6799999999999999e-012	1	1	\N	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	4.5	\N
13	{"format": "CCP4", "sizeKb": 364501, "type": "map", "file": "emd_0010.map.gz"}	Image stored as Reals	450	450	450	0	0	0	449	449	449	450	450	450	{"units": "A", "value": 495.0}	{"units": "A", "value": 495.0}	{"units": "A", "value": 495.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	Z	Y	X	-13.500279000000001	32.492694999999998	3.6600000000000001e-013	1	1	\N	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	6	\N
14	{"format": "CCP4", "sizeKb": 256001, "type": "map", "file": "emd_0011.map.gz"}	Image stored as Reals	400	400	400	0	0	0	399	399	399	400	400	400	{"units": "A", "value": 421.6}	{"units": "A", "value": 421.6}	{"units": "A", "value": 421.6}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.067202910000000005	0.15161994000000001	0.00038904152	0.0052055469999999996	1	\N	{"units": "A", "value": 1.054}	{"units": "A", "value": 1.054}	{"units": "A", "value": 1.054}	0.021000000000000001	\N
15	{"format": "CCP4", "sizeKb": 44958, "type": "map", "file": "emd_0012.map.gz"}	Image stored as Reals	224	224	224	0	0	0	223	223	223	224	224	224	{"units": "A", "value": 239.68001}	{"units": "A", "value": 239.68001}	{"units": "A", "value": 239.68001}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.075805224000000004	0.18155326999999999	0.0010824882	0.0091030250000000007	1	\N	{"units": "A", "value": 1.07}	{"units": "A", "value": 1.07}	{"units": "A", "value": 1.07}	0.049500000000000002	MDA5-dsRNA helical reconstruction in the presence of 2.5 mM AMPPNP
16	{"format": "CCP4", "sizeKb": 96101, "type": "map", "file": "emd_0013.map.gz"}	Image stored as Reals	310	310	250	-155	-155	-125	154	154	124	310	310	250	{"units": "A", "value": 310.62}	{"units": "A", "value": 310.62}	{"units": "A", "value": 250.5}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-3.8562029999999998	7.4975509999999996	0.019812632	0.63959557	1	\N	{"units": "A", "value": 1.002}	{"units": "A", "value": 1.002}	{"units": "A", "value": 1.002}	2	\N
17	{"format": "CCP4", "sizeKb": 108001, "type": "map", "file": "emd_0014.map.gz"}	Image stored as Reals	300	300	300	0	0	0	299	299	299	300	300	300	{"units": "A", "value": 317.99997}	{"units": "A", "value": 317.99997}	{"units": "A", "value": 317.99997}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.081275730000000004	0.13223831	0.00026691699999999998	0.0030406867000000001	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.039800000000000002	Beta-2-microglobulin amyloid fibril - two protofilamets
18	{"format": "CCP4", "sizeKb": 67109, "type": "map", "file": "emd_0015.map.gz"}	Image stored as Reals	256	256	256	0	0	0	255	255	255	256	256	256	{"units": "A", "value": 473.6}	{"units": "A", "value": 473.6}	{"units": "A", "value": 473.6}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.11295794000000001	0.27922671999999998	0.0012513183000000001	0.013090289999999999	1	\N	{"units": "A", "value": 1.85}	{"units": "A", "value": 1.85}	{"units": "A", "value": 1.85}	0.079699999999999993	\N
19	{"format": "CCP4", "sizeKb": 67109, "type": "map", "file": "emd_0001.map.gz"}	Image stored as Reals	256	256	256	0	0	0	255	255	255	256	256	256	{"units": "A", "value": 1.0}	{"units": "A", "value": 1.0}	{"units": "A", "value": 1.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.10788548000000001	0.22688343999999999	1.2941896e-005	0.0077065425	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.017999999999999999	\N
20	{"format": "CCP4", "sizeKb": 67109, "type": "map", "file": "emd_0002.map.gz"}	Image stored as Reals	256	256	256	0	0	0	255	255	255	256	256	256	{"units": "A", "value": 276.48}	{"units": "A", "value": 276.48}	{"units": "A", "value": 276.48}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.081265143999999997	0.19144611	0.00053819930000000005	0.0061218822000000004	1	\N	{"units": "A", "value": 1.08}	{"units": "A", "value": 1.08}	{"units": "A", "value": 1.08}	0.016	\N
21	{"format": "CCP4", "sizeKb": 202613, "type": "map", "file": "emd_0004.map.gz"}	Image stored as Reals	370	370	370	0	0	0	369	369	369	370	370	370	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-1.2225216999999999	2.4226390000000002	0.0082020170000000007	0.080526760000000003	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.33000000000000002	Modelling map
22	{"format": "CCP4", "sizeKb": 202613, "type": "map", "file": "emd_0005.map.gz"}	Image stored as Reals	370	370	370	0	0	0	369	369	369	370	370	370	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-1.0410119	2.3215709000000002	0.0087016790000000004	0.076787889999999998	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.33000000000000002	Saccharomyces cerevisiae Class1 (III2IV)
23	{"format": "CCP4", "sizeKb": 202613, "type": "map", "file": "emd_0006.map.gz"}	Image stored as Reals	370	370	370	0	0	0	369	369	369	370	370	370	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "A", "value": 392.19998}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.57250069999999997	1.5357721	0.0071697747000000001	0.061136416999999998	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.22	Saccharomyces cerevisiae supercomplex class2 (III2IV2)
24	{"format": "CCP4", "sizeKb": 32001, "type": "map", "file": "emd_0007.map.gz"}	Image stored as Reals	200	200	200	0	0	0	199	199	199	200	200	200	{"units": "A", "value": 252.0}	{"units": "A", "value": 252.0}	{"units": "A", "value": 252.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.026832027000000001	0.10347184500000001	0.0019379637	0.0087864729999999995	1	\N	{"units": "A", "value": 1.26}	{"units": "A", "value": 1.26}	{"units": "A", "value": 1.26}	0.042999999999999997	\N
25	{"format": "CCP4", "sizeKb": 364501, "type": "map", "file": "emd_0008.map.gz"}	Image stored as Reals	450	450	450	0	0	0	449	449	449	450	450	450	{"units": "A", "value": 495.0}	{"units": "A", "value": 495.0}	{"units": "A", "value": 495.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	Z	Y	X	-11.851832	30.160060000000001	4.4800000000000003e-012	1	1	\N	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	6.5	\N
26	{"format": "CCP4", "sizeKb": 32001, "type": "map", "file": "emd_0009.map.gz"}	Image stored as Reals	200	200	200	0	0	0	199	199	199	200	200	200	{"units": "A", "value": 220.0}	{"units": "A", "value": 220.0}	{"units": "A", "value": 220.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	Z	Y	X	-9.0799559999999992	16.808064999999999	-2.6799999999999999e-012	1	1	\N	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	4.5	\N
27	{"format": "CCP4", "sizeKb": 364501, "type": "map", "file": "emd_0010.map.gz"}	Image stored as Reals	450	450	450	0	0	0	449	449	449	450	450	450	{"units": "A", "value": 495.0}	{"units": "A", "value": 495.0}	{"units": "A", "value": 495.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	Z	Y	X	-13.500279000000001	32.492694999999998	3.6600000000000001e-013	1	1	\N	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	{"units": "A", "value": 1.1}	6	\N
28	{"format": "CCP4", "sizeKb": 256001, "type": "map", "file": "emd_0011.map.gz"}	Image stored as Reals	400	400	400	0	0	0	399	399	399	400	400	400	{"units": "A", "value": 421.6}	{"units": "A", "value": 421.6}	{"units": "A", "value": 421.6}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.067202910000000005	0.15161994000000001	0.00038904152	0.0052055469999999996	1	\N	{"units": "A", "value": 1.054}	{"units": "A", "value": 1.054}	{"units": "A", "value": 1.054}	0.021000000000000001	\N
29	{"format": "CCP4", "sizeKb": 44958, "type": "map", "file": "emd_0012.map.gz"}	Image stored as Reals	224	224	224	0	0	0	223	223	223	224	224	224	{"units": "A", "value": 239.68001}	{"units": "A", "value": 239.68001}	{"units": "A", "value": 239.68001}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.075805224000000004	0.18155326999999999	0.0010824882	0.0091030250000000007	1	\N	{"units": "A", "value": 1.07}	{"units": "A", "value": 1.07}	{"units": "A", "value": 1.07}	0.049500000000000002	MDA5-dsRNA helical reconstruction in the presence of 2.5 mM AMPPNP
30	{"format": "CCP4", "sizeKb": 96101, "type": "map", "file": "emd_0013.map.gz"}	Image stored as Reals	310	310	250	-155	-155	-125	154	154	124	310	310	250	{"units": "A", "value": 310.62}	{"units": "A", "value": 310.62}	{"units": "A", "value": 250.5}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-3.8562029999999998	7.4975509999999996	0.019812632	0.63959557	1	\N	{"units": "A", "value": 1.002}	{"units": "A", "value": 1.002}	{"units": "A", "value": 1.002}	2	\N
31	{"format": "CCP4", "sizeKb": 108001, "type": "map", "file": "emd_0014.map.gz"}	Image stored as Reals	300	300	300	0	0	0	299	299	299	300	300	300	{"units": "A", "value": 317.99997}	{"units": "A", "value": 317.99997}	{"units": "A", "value": 317.99997}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.081275730000000004	0.13223831	0.00026691699999999998	0.0030406867000000001	1	\N	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	{"units": "A", "value": 1.06}	0.039800000000000002	Beta-2-microglobulin amyloid fibril - two protofilamets
32	{"format": "CCP4", "sizeKb": 67109, "type": "map", "file": "emd_0015.map.gz"}	Image stored as Reals	256	256	256	0	0	0	255	255	255	256	256	256	{"units": "A", "value": 473.6}	{"units": "A", "value": 473.6}	{"units": "A", "value": 473.6}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.11295794000000001	0.27922671999999998	0.0012513183000000001	0.013090289999999999	1	\N	{"units": "A", "value": 1.85}	{"units": "A", "value": 1.85}	{"units": "A", "value": 1.85}	0.079699999999999993	\N
33	{"format": "CCP4", "sizeKb": 78733, "type": "map", "file": "emd_0016.map.gz"}	Image stored as Reals	270	270	270	0	0	0	269	269	269	270	270	270	{"units": "A", "value": 499.5}	{"units": "A", "value": 499.5}	{"units": "A", "value": 499.5}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.086722770000000005	0.30624499999999999	0.0013167954000000001	0.013812984	1	\N	{"units": "A", "value": 1.85}	{"units": "A", "value": 1.85}	{"units": "A", "value": 1.85}	0.071999999999999995	\N
34	{"format": "CCP4", "sizeKb": 78733, "type": "map", "file": "emd_0017.map.gz"}	Image stored as Reals	270	270	270	0	0	0	269	269	269	270	270	270	{"units": "A", "value": 499.5}	{"units": "A", "value": 499.5}	{"units": "A", "value": 499.5}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.082781859999999999	0.27012946999999998	0.0012991953000000001	0.013124164000000001	1	\N	{"units": "A", "value": 1.85}	{"units": "A", "value": 1.85}	{"units": "A", "value": 1.85}	0.054899999999999997	\N
35	{"format": "CCP4", "sizeKb": 78733, "type": "map", "file": "emd_0018.map.gz"}	Image stored as Reals	270	270	270	0	0	0	269	269	269	270	270	270	{"units": "A", "value": 499.5}	{"units": "A", "value": 499.5}	{"units": "A", "value": 499.5}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.026710883000000001	0.084918049999999995	9.24297e-005	0.0057740485999999997	1	\N	{"units": "A", "value": 1.85}	{"units": "A", "value": 1.85}	{"units": "A", "value": 1.85}	0.0043	\N
36	{"format": "CCP4", "sizeKb": 67109, "type": "map", "file": "emd_0020.map.gz"}	Image stored as Reals	256	256	256	0	0	0	255	255	255	256	256	256	{"units": "A", "value": 220.16}	{"units": "A", "value": 220.16}	{"units": "A", "value": 220.16}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	{"units": "degrees", "value": 90.0}	X	Y	Z	-0.92154700000000001	1.6588855	0.0027781222000000001	0.061088625000000001	1	\N	{"units": "A", "value": 0.86}	{"units": "A", "value": 0.86}	{"units": "A", "value": 0.86}	0.23000000000000001	Cytochrome c nitrite reductase from the bacterium Thioalkalivibrio nitratireducens (TvNiR).
\.


--
-- TOC entry 2932 (class 0 OID 16532)
-- Dependencies: 200
-- Data for Name: emd_entry; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.emd_entry (id, full_name, acronym, volume, resolution, image_url, xml_url, map_url, map_information_id, png_img_3d, gif_img_3d) FROM stdin;
1	Cryo-EM structure of bacterial RNA polymerase-sigma54 holoenzyme transcription open complex	Cryo-EM structure of bacterial RNA polymerase-sigma54 holoenzyme transcription open complex	70	3.3999999999999999	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0001/images/emd_0001.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0001/header/emd-0001.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0001/map/emd_0001.map.gz	19	assets/img/test-3d.png	assets/img/test-3d.gif
2	Cryo-EM structure of bacterial RNA polymerase-sigma54 holoenzyme intermediate partially loaded complex	Cryo-EM structure of bacterial RNA polymerase-sigma54 holoenzyme intermediate partially loaded complex	70	4.0999999999999996	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0002/images/emd_0002.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0002/header/emd-0002.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0002/map/emd_0002.map.gz	20	assets/img/test-3d.png	assets/img/test-3d.gif
4	Saccharomyces cerevisiae Supercomplex (III2IV)	Saccharomyces cerevisiae Supercomplex (III2IV)	26	3.23	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0004/images/emd_0004.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0004/header/emd-0004.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0004/map/emd_0004.map.gz	21	assets/img/test-3d.png	assets/img/test-3d.gif
5	Saccharomyces cerevisiae Supercomplex (III2IV)	Saccharomyces cerevisiae Supercomplex (III2IV)	26	3.3399999999999999	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0005/images/emd_0005.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0005/header/emd-0005.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0005/map/emd_0005.map.gz	22	assets/img/test-3d.png	assets/img/test-3d.gif
6	Saccharomyces cerevisiae Supercomplex	Saccharomyces cerevisiae Supercomplex	26	3.5	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0006/images/emd_0006.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0006/header/emd-0006.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0006/map/emd_0006.map.gz	23	assets/img/test-3d.png	assets/img/test-3d.gif
7	MlaBDEF	MlaBDEF	8	8.6999999999999993	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0007/images/emd_0007.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0007/header/emd-0007.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0007/map/emd_0007.map.gz	24	assets/img/test-3d.png	assets/img/test-3d.gif
8	Complex of TssK, TssF and TssG, components of the type VI secretion system baseplate	Complex of TssK, TssF and TssG, components of the type VI secretion system baseplate	3	4.2999999999999998	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0008/images/emd_0008.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0008/header/emd-0008.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0008/map/emd_0008.map.gz	25	assets/img/test-3d.png	assets/img/test-3d.gif
9	Type VI secretion system baseplate complex	Type VI secretion system baseplate complex	3	4.7000000000000002	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0009/images/emd_0009.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0009/header/emd-0009.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0009/map/emd_0009.map.gz	26	assets/img/test-3d.png	assets/img/test-3d.gif
10	Complex of TssK, TssF and TssG, components of the type VI secretion system baseplate	Complex of TssK, TssF and TssG, components of the type VI secretion system baseplate	3	4.2999999999999998	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0010/images/emd_0010.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0010/header/emd-0010.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0010/map/emd_0010.map.gz	27	assets/img/test-3d.png	assets/img/test-3d.gif
11	Fatty Acid Synthase - I	Fatty Acid Synthase - I	9	3.2999999999999998	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0011/images/emd_0011.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0011/header/emd-0011.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0011/map/emd_0011.map.gz	28	assets/img/test-3d.png	assets/img/test-3d.gif
12	MDA5-dsRNA helical filament in complex with AMPPNP	MDA5-dsRNA helical filament in complex with AMPPNP	72	4.0599999999999996	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0012/images/emd_0012.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0012/header/emd-0012.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0012/map/emd_0012.map.gz	29	assets/img/test-3d.png	assets/img/test-3d.gif
13	Complex of BCL10 CARD and MALT1 DEATH DOMAIN	Complex of BCL10 CARD and MALT1 DEATH DOMAIN	9	4.9000000000000004	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0013/images/emd_0013.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0013/header/emd-0013.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0013/map/emd_0013.map.gz	30	assets/img/test-3d.png	assets/img/test-3d.gif
14	two protofilament beta-2-microglobulin amyloid fibril	two protofilament beta-2-microglobulin amyloid fibril	9	3.9750000000000001	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0014/images/emd_0014.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0014/header/emd-0014.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0014/map/emd_0014.map.gz	31	assets/img/test-3d.png	assets/img/test-3d.gif
15	GroEL in complex with actin	GroEL in complex with actin	174	10.300000000000001	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0015/images/emd_0015.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0015/header/emd-0015.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0015/map/emd_0015.map.gz	32	assets/img/test-3d.png	assets/img/test-3d.gif
16	Bovine TRiC in nucleotide free state	Bovine TRiC in nucleotide free state	174	8.8000000000000007	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0016/images/emd_0016.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0016/header/emd-0016.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0016/map/emd_0016.map.gz	33	assets/img/test-3d.png	assets/img/test-3d.gif
17	Bovine TRiC in complex with actin	Bovine TRiC in complex with actin	174	9.0999999999999996	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0017/images/emd_0017.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0017/header/emd-0017.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0017/map/emd_0017.map.gz	34	assets/img/test-3d.png	assets/img/test-3d.gif
18	Bovine TRiC in complex with actin and CCT1 antibody	Bovine TRiC in complex with actin and CCT1 antibody	174	16.699999999999999	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0018/images/emd_0018.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0018/header/emd-0018.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0018/map/emd_0018.map.gz	35	assets/img/test-3d.png	assets/img/test-3d.gif
20	Cytochrome c nitrite reductase from the bacterium Thioalkalivibrio nitratireducens (TvNiR)	Cytochrome c nitrite reductase from the bacterium Thioalkalivibrio nitratireducens (TvNiR)	10	2.5600000000000001	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0020/images/emd_0020.png	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0020/header/emd-0020.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0020/map/emd_0020.map.gz	36	assets/img/test-3d.png	assets/img/test-3d.gif
\.


--
-- TOC entry 2931 (class 0 OID 16521)
-- Dependencies: 199
-- Data for Name: type_descriptor; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.type_descriptor (id, name, description) FROM stdin;
0	test	test
\.


--
-- TOC entry 2933 (class 0 OID 16545)
-- Dependencies: 201
-- Data for Name: descriptor; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.descriptor (emd_entry_id, type_descriptor_id, numbers) FROM stdin;
1	0	{"descriptors":[1,2,3,4,5]}
2	0	[1,2,3,4,5]
4	0	[6,7,8,9]
10	0	[45.6,34.7,96.3,65.8,9.12]
7	0	[12,60,79,69,14,5,36,17,35,59,65,70,58,94,3,23,37,9,34,77,39,73,81,87,92,68,100,96,89,86,62,16,38,8,78,33,66,52,83,53,84,44,99,91,76,54,82,30,80,25,55,26,19,97,20,22,11,7,57,61,1,15,67,50,24,18,43,49]
14	0	[96,89,86,62,16,38,8,78,33,66,52,83,53,84,44,99,91,76,54,82,30,80,25,55,26,19,12,60,79,69,14,5,36,17,35,59,65,70,58,94,3,23,37,9,34,77,39,73,81,87,92,68,100,97,20,22,11,7,57,61,1,15,67,50,24,18,43,49]
\.


--
-- TOC entry 2939 (class 0 OID 16617)
-- Dependencies: 207
-- Data for Name: pdb_entry; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.pdb_entry (id, pdb) FROM stdin;
\.


--
-- TOC entry 2940 (class 0 OID 16626)
-- Dependencies: 208
-- Data for Name: pdb_entry_x_emd_entry; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.pdb_entry_x_emd_entry (pdb_entry_id, emd_entry_id) FROM stdin;
\.


--
-- TOC entry 2949 (class 0 OID 17035)
-- Dependencies: 217
-- Data for Name: search_history; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.search_history (id, date_time, ip, user_id, emd_entry_id, name_file, contour_level, representation_id, volume_filter_id, resolution_filter_min, resolution_filter_max) FROM stdin;
1	2019-05-24 00:56:34.932	::1	\N	14	-	\N	2	0	3	9
2	2019-05-24 00:57:11.65	::1	\N	\N	ejemplo.map	3.1400000000000001	4	0	2	78
3	2019-05-24 01:05:45.196	::1	\N	4	-	\N	2	0	1	15
\.


--
-- TOC entry 2941 (class 0 OID 16641)
-- Dependencies: 209
-- Data for Name: time_stamp; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.time_stamp (emd_entry_id, modification, map_file, xml_file, image_file) FROM stdin;
1	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0001/map/emd_0001.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0001/header/emd-0001.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0001/images/emd_0001.png
2	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0002/map/emd_0002.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0002/header/emd-0002.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0002/images/emd_0002.png
4	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0004/map/emd_0004.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0004/header/emd-0004.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0004/images/emd_0004.png
5	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0005/map/emd_0005.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0005/header/emd-0005.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0005/images/emd_0005.png
6	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0006/map/emd_0006.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0006/header/emd-0006.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0006/images/emd_0006.png
7	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0007/map/emd_0007.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0007/header/emd-0007.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0007/images/emd_0007.png
8	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0008/map/emd_0008.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0008/header/emd-0008.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0008/images/emd_0008.png
9	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0009/map/emd_0009.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0009/header/emd-0009.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0009/images/emd_0009.png
10	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0010/map/emd_0010.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0010/header/emd-0010.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0010/images/emd_0010.png
11	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0011/map/emd_0011.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0011/header/emd-0011.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0011/images/emd_0011.png
12	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0012/map/emd_0012.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0012/header/emd-0012.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0012/images/emd_0012.png
13	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0013/map/emd_0013.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0013/header/emd-0013.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0013/images/emd_0013.png
14	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0014/map/emd_0014.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0014/header/emd-0014.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0014/images/emd_0014.png
15	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0015/map/emd_0015.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0015/header/emd-0015.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0015/images/emd_0015.png
16	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0016/map/emd_0016.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0016/header/emd-0016.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0016/images/emd_0016.png
17	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0017/map/emd_0017.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0017/header/emd-0017.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0017/images/emd_0017.png
18	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0018/map/emd_0018.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0018/header/emd-0018.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0018/images/emd_0018.png
20	2019-04-07	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0020/map/emd_0020.map.gz	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0020/header/emd-0020.xml	http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0020/images/emd_0020.png
\.


--
-- TOC entry 2942 (class 0 OID 16654)
-- Dependencies: 210
-- Data for Name: update; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.update (last_update) FROM stdin;
2019-04-07
\.


--
-- TOC entry 2955 (class 0 OID 0)
-- Dependencies: 214
-- Name: benchmark_history_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.benchmark_history_id_seq', 2, true);


--
-- TOC entry 2956 (class 0 OID 0)
-- Dependencies: 196
-- Name: map_information_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.map_information_id_seq', 36, true);


--
-- TOC entry 2957 (class 0 OID 0)
-- Dependencies: 206
-- Name: pdb_entry_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.pdb_entry_id_seq', 1, false);


--
-- TOC entry 2958 (class 0 OID 0)
-- Dependencies: 202
-- Name: representation_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.representation_id_seq', 1, false);


--
-- TOC entry 2959 (class 0 OID 0)
-- Dependencies: 216
-- Name: search_history_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.search_history_id_seq', 3, true);


--
-- TOC entry 2960 (class 0 OID 0)
-- Dependencies: 198
-- Name: type_descriptor_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.type_descriptor_id_seq', 1, false);


--
-- TOC entry 2961 (class 0 OID 0)
-- Dependencies: 211
-- Name: user_role_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.user_role_id_seq', 1, false);


--
-- TOC entry 2962 (class 0 OID 0)
-- Dependencies: 204
-- Name: volume_filter_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.volume_filter_id_seq', 1, false);


-- Completed on 2019-05-23 19:10:43

--
-- PostgreSQL database dump complete
--

