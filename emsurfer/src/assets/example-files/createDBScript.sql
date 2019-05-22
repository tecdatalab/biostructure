--
-- PostgreSQL database dump
--

-- Dumped from database version 10.7
-- Dumped by pg_dump version 10.7

-- Started on 2019-05-21 14:18:51

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
-- TOC entry 1 (class 3079 OID 12924)
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- TOC entry 2975 (class 0 OID 0)
-- Dependencies: 1
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


SET default_tablespace = '';

SET default_with_oids = false;

--
-- TOC entry 212 (class 1259 OID 16704)
-- Name: benchmark_history; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.benchmark_history (
    id integer NOT NULL,
    date_time timestamp without time zone NOT NULL,
    ip text NOT NULL,
    user_id integer,
    representation_id integer NOT NULL,
    volume_filter_id integer NOT NULL,
    top_results integer NOT NULL,
    emd_list json NOT NULL
);


ALTER TABLE public.benchmark_history OWNER TO postgres;

--
-- TOC entry 2976 (class 0 OID 0)
-- Dependencies: 212
-- Name: TABLE benchmark_history; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.benchmark_history IS 'Stores the significant data of the benchmark queries.';


--
-- TOC entry 211 (class 1259 OID 16702)
-- Name: benchmark_history_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.benchmark_history_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.benchmark_history_id_seq OWNER TO postgres;

--
-- TOC entry 2977 (class 0 OID 0)
-- Dependencies: 211
-- Name: benchmark_history_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.benchmark_history_id_seq OWNED BY public.benchmark_history.id;


--
-- TOC entry 221 (class 1259 OID 16852)
-- Name: benchmarks_history; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.benchmarks_history (
    id integer NOT NULL,
    date_time timestamp without time zone NOT NULL,
    ip text NOT NULL,
    user_id text NOT NULL,
    representation_id integer NOT NULL,
    volume_filter_id integer NOT NULL,
    top_results integer NOT NULL,
    emd_list json NOT NULL
);


ALTER TABLE public.benchmarks_history OWNER TO postgres;

--
-- TOC entry 220 (class 1259 OID 16850)
-- Name: benchmarks_history_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.benchmarks_history_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.benchmarks_history_id_seq OWNER TO postgres;

--
-- TOC entry 2978 (class 0 OID 0)
-- Dependencies: 220
-- Name: benchmarks_history_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.benchmarks_history_id_seq OWNED BY public.benchmarks_history.id;


--
-- TOC entry 201 (class 1259 OID 16545)
-- Name: descriptor; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.descriptor (
    emd_entry_id integer NOT NULL,
    type_descriptor_id integer NOT NULL,
    numbers json NOT NULL
);


ALTER TABLE public.descriptor OWNER TO postgres;

--
-- TOC entry 2979 (class 0 OID 0)
-- Dependencies: 201
-- Name: TABLE descriptor; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.descriptor IS 'It relates the EMD with the descriptors, in a one-to-one relationship.';


--
-- TOC entry 200 (class 1259 OID 16532)
-- Name: emd_entry; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.emd_entry (
    id integer NOT NULL,
    full_name text NOT NULL,
    acronym text NOT NULL,
    volume double precision NOT NULL,
    resolution double precision NOT NULL,
    image_url text,
    xml_url text NOT NULL,
    map_url text NOT NULL,
    map_information_id integer NOT NULL,
    png_img_3d text,
    gif_img_3d text
);


ALTER TABLE public.emd_entry OWNER TO postgres;

--
-- TOC entry 2980 (class 0 OID 0)
-- Dependencies: 200
-- Name: TABLE emd_entry; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.emd_entry IS 'Stores existing EMDs.';


--
-- TOC entry 197 (class 1259 OID 16510)
-- Name: map_information; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.map_information (
    id integer NOT NULL,
    file_information json,
    data_type text,
    num_columns integer,
    num_rows integer,
    num_sections integer,
    origin_col integer,
    origin_row integer,
    origin_sec integer,
    limit_col integer,
    limit_row integer,
    limit_sec integer,
    spacing_col integer,
    spacing_row integer,
    spacing_sec integer,
    cell_a json,
    cell_b json,
    cell_c json,
    cell_alpha json,
    cell_beta json,
    cell_gamma json,
    axis_order_fast character(1),
    axis_order_medium character(1),
    axis_order_slow character(1),
    minimum double precision,
    maximum double precision,
    average double precision,
    std double precision,
    space_group_number integer,
    details text,
    pixel_x json,
    pixel_y json,
    pixel_z json,
    countour_level double precision,
    annotation_details text
);


ALTER TABLE public.map_information OWNER TO postgres;

--
-- TOC entry 2981 (class 0 OID 0)
-- Dependencies: 197
-- Name: TABLE map_information; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.map_information IS 'Stores the map information of an EMD.';


--
-- TOC entry 196 (class 1259 OID 16508)
-- Name: map_information_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.map_information_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.map_information_id_seq OWNER TO postgres;

--
-- TOC entry 2982 (class 0 OID 0)
-- Dependencies: 196
-- Name: map_information_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.map_information_id_seq OWNED BY public.map_information.id;


--
-- TOC entry 207 (class 1259 OID 16617)
-- Name: pdb_entry; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.pdb_entry (
    id integer NOT NULL,
    pdb text NOT NULL
);


ALTER TABLE public.pdb_entry OWNER TO postgres;

--
-- TOC entry 2983 (class 0 OID 0)
-- Dependencies: 207
-- Name: TABLE pdb_entry; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.pdb_entry IS 'Stores the different PDBs.';


--
-- TOC entry 206 (class 1259 OID 16615)
-- Name: pdb_entry_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.pdb_entry_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.pdb_entry_id_seq OWNER TO postgres;

--
-- TOC entry 2984 (class 0 OID 0)
-- Dependencies: 206
-- Name: pdb_entry_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.pdb_entry_id_seq OWNED BY public.pdb_entry.id;


--
-- TOC entry 208 (class 1259 OID 16626)
-- Name: pdb_entry_x_emd_entry; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.pdb_entry_x_emd_entry (
    pdb_entry_id integer NOT NULL,
    emd_entry_id integer NOT NULL
);


ALTER TABLE public.pdb_entry_x_emd_entry OWNER TO postgres;

--
-- TOC entry 2985 (class 0 OID 0)
-- Dependencies: 208
-- Name: TABLE pdb_entry_x_emd_entry; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.pdb_entry_x_emd_entry IS 'Associate PDBs with EMDs in a many-to-many relationship.';


--
-- TOC entry 203 (class 1259 OID 16565)
-- Name: representation; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.representation (
    id integer NOT NULL,
    name text NOT NULL
);


ALTER TABLE public.representation OWNER TO postgres;

--
-- TOC entry 2986 (class 0 OID 0)
-- Dependencies: 203
-- Name: TABLE representation; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.representation IS 'Stores the types of representations of the contour shape.';


--
-- TOC entry 202 (class 1259 OID 16563)
-- Name: representation_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.representation_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.representation_id_seq OWNER TO postgres;

--
-- TOC entry 2987 (class 0 OID 0)
-- Dependencies: 202
-- Name: representation_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.representation_id_seq OWNED BY public.representation.id;


--
-- TOC entry 223 (class 1259 OID 16878)
-- Name: search_history; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.search_history (
    id integer NOT NULL,
    date_time timestamp without time zone NOT NULL,
    ip text NOT NULL,
    user_id text NOT NULL,
    emd_entry_id integer NOT NULL,
    name_file text,
    counter_level double precision,
    representation_id integer NOT NULL,
    volume_filter_id integer NOT NULL,
    resolution_filter_min double precision,
    resolution_filter_max double precision
);


ALTER TABLE public.search_history OWNER TO postgres;

--
-- TOC entry 222 (class 1259 OID 16876)
-- Name: search_history_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.search_history_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.search_history_id_seq OWNER TO postgres;

--
-- TOC entry 2988 (class 0 OID 0)
-- Dependencies: 222
-- Name: search_history_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.search_history_id_seq OWNED BY public.search_history.id;


--
-- TOC entry 209 (class 1259 OID 16641)
-- Name: time_stamp; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.time_stamp (
    emd_entry_id integer NOT NULL,
    modification date NOT NULL,
    map_file text NOT NULL,
    xml_file text NOT NULL,
    image_file text NOT NULL
);


ALTER TABLE public.time_stamp OWNER TO postgres;

--
-- TOC entry 2989 (class 0 OID 0)
-- Dependencies: 209
-- Name: TABLE time_stamp; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.time_stamp IS 'Table to know the date of the last EMD update of the database.';


--
-- TOC entry 199 (class 1259 OID 16521)
-- Name: type_descriptor; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.type_descriptor (
    id integer NOT NULL,
    name text NOT NULL,
    description text NOT NULL
);


ALTER TABLE public.type_descriptor OWNER TO postgres;

--
-- TOC entry 2990 (class 0 OID 0)
-- Dependencies: 199
-- Name: TABLE type_descriptor; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.type_descriptor IS 'Stores the different types of existing descriptors.';


--
-- TOC entry 198 (class 1259 OID 16519)
-- Name: type_descriptor_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.type_descriptor_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.type_descriptor_id_seq OWNER TO postgres;

--
-- TOC entry 2991 (class 0 OID 0)
-- Dependencies: 198
-- Name: type_descriptor_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.type_descriptor_id_seq OWNED BY public.type_descriptor.id;


--
-- TOC entry 210 (class 1259 OID 16654)
-- Name: update; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.update (
    last_update date NOT NULL
);


ALTER TABLE public.update OWNER TO postgres;

--
-- TOC entry 2992 (class 0 OID 0)
-- Dependencies: 210
-- Name: TABLE update; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.update IS 'Table to know the last day in which the updating procedures were executed.';


--
-- TOC entry 219 (class 1259 OID 16808)
-- Name: user; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public."user" (
    id text NOT NULL,
    name text NOT NULL,
    email text NOT NULL,
    role integer DEFAULT 1 NOT NULL
);


ALTER TABLE public."user" OWNER TO postgres;

--
-- TOC entry 2993 (class 0 OID 0)
-- Dependencies: 219
-- Name: TABLE "user"; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public."user" IS 'Stores the user information';


--
-- TOC entry 218 (class 1259 OID 16780)
-- Name: user2; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.user2 (
    id text NOT NULL,
    name text NOT NULL,
    email text NOT NULL,
    role integer DEFAULT 0 NOT NULL
);


ALTER TABLE public.user2 OWNER TO postgres;

--
-- TOC entry 216 (class 1259 OID 16757)
-- Name: user_role; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.user_role (
    id integer NOT NULL,
    role text NOT NULL
);


ALTER TABLE public.user_role OWNER TO postgres;

--
-- TOC entry 215 (class 1259 OID 16755)
-- Name: user_role_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.user_role_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.user_role_id_seq OWNER TO postgres;

--
-- TOC entry 2994 (class 0 OID 0)
-- Dependencies: 215
-- Name: user_role_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.user_role_id_seq OWNED BY public.user_role.id;


--
-- TOC entry 214 (class 1259 OID 16725)
-- Name: users_roles; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.users_roles (
    id integer NOT NULL,
    rol text NOT NULL
);


ALTER TABLE public.users_roles OWNER TO postgres;

--
-- TOC entry 2995 (class 0 OID 0)
-- Dependencies: 214
-- Name: TABLE users_roles; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.users_roles IS 'Stores the different types of user roles ';


--
-- TOC entry 213 (class 1259 OID 16723)
-- Name: users_roles_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.users_roles_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.users_roles_id_seq OWNER TO postgres;

--
-- TOC entry 2996 (class 0 OID 0)
-- Dependencies: 213
-- Name: users_roles_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.users_roles_id_seq OWNED BY public.users_roles.id;


--
-- TOC entry 217 (class 1259 OID 16766)
-- Name: userss; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.userss (
    id text NOT NULL,
    name text NOT NULL,
    email text NOT NULL,
    role integer DEFAULT 0 NOT NULL
);


ALTER TABLE public.userss OWNER TO postgres;

--
-- TOC entry 205 (class 1259 OID 16578)
-- Name: volume_filter; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.volume_filter (
    id integer NOT NULL,
    name text NOT NULL
);


ALTER TABLE public.volume_filter OWNER TO postgres;

--
-- TOC entry 2997 (class 0 OID 0)
-- Dependencies: 205
-- Name: TABLE volume_filter; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.volume_filter IS 'Stores the states that the volume filter can be (off, on).';


--
-- TOC entry 204 (class 1259 OID 16576)
-- Name: volume_filter_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.volume_filter_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.volume_filter_id_seq OWNER TO postgres;

--
-- TOC entry 2998 (class 0 OID 0)
-- Dependencies: 204
-- Name: volume_filter_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.volume_filter_id_seq OWNED BY public.volume_filter.id;


--
-- TOC entry 2777 (class 2604 OID 16707)
-- Name: benchmark_history id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.benchmark_history ALTER COLUMN id SET DEFAULT nextval('public.benchmark_history_id_seq'::regclass);


--
-- TOC entry 2783 (class 2604 OID 16855)
-- Name: benchmarks_history id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.benchmarks_history ALTER COLUMN id SET DEFAULT nextval('public.benchmarks_history_id_seq'::regclass);


--
-- TOC entry 2772 (class 2604 OID 16513)
-- Name: map_information id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.map_information ALTER COLUMN id SET DEFAULT nextval('public.map_information_id_seq'::regclass);


--
-- TOC entry 2776 (class 2604 OID 16620)
-- Name: pdb_entry id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pdb_entry ALTER COLUMN id SET DEFAULT nextval('public.pdb_entry_id_seq'::regclass);


--
-- TOC entry 2774 (class 2604 OID 16568)
-- Name: representation id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.representation ALTER COLUMN id SET DEFAULT nextval('public.representation_id_seq'::regclass);


--
-- TOC entry 2784 (class 2604 OID 16881)
-- Name: search_history id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.search_history ALTER COLUMN id SET DEFAULT nextval('public.search_history_id_seq'::regclass);


--
-- TOC entry 2773 (class 2604 OID 16524)
-- Name: type_descriptor id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.type_descriptor ALTER COLUMN id SET DEFAULT nextval('public.type_descriptor_id_seq'::regclass);


--
-- TOC entry 2779 (class 2604 OID 16760)
-- Name: user_role id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.user_role ALTER COLUMN id SET DEFAULT nextval('public.user_role_id_seq'::regclass);


--
-- TOC entry 2778 (class 2604 OID 16728)
-- Name: users_roles id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.users_roles ALTER COLUMN id SET DEFAULT nextval('public.users_roles_id_seq'::regclass);


--
-- TOC entry 2775 (class 2604 OID 16581)
-- Name: volume_filter id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.volume_filter ALTER COLUMN id SET DEFAULT nextval('public.volume_filter_id_seq'::regclass);


--
-- TOC entry 2812 (class 2606 OID 16712)
-- Name: benchmark_history benchmark_history_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.benchmark_history
    ADD CONSTRAINT benchmark_history_pkey PRIMARY KEY (id);


--
-- TOC entry 2826 (class 2606 OID 16860)
-- Name: benchmarks_history benchmarks_history_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.benchmarks_history
    ADD CONSTRAINT benchmarks_history_pkey PRIMARY KEY (id);


--
-- TOC entry 2794 (class 2606 OID 16552)
-- Name: descriptor descriptor_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.descriptor
    ADD CONSTRAINT descriptor_pkey PRIMARY KEY (emd_entry_id, type_descriptor_id);


--
-- TOC entry 2792 (class 2606 OID 16539)
-- Name: emd_entry emd_entry_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.emd_entry
    ADD CONSTRAINT emd_entry_pkey PRIMARY KEY (id);


--
-- TOC entry 2786 (class 2606 OID 16518)
-- Name: map_information map_information_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.map_information
    ADD CONSTRAINT map_information_pkey PRIMARY KEY (id);


--
-- TOC entry 2804 (class 2606 OID 16625)
-- Name: pdb_entry pdb_entry_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pdb_entry
    ADD CONSTRAINT pdb_entry_pkey PRIMARY KEY (id);


--
-- TOC entry 2806 (class 2606 OID 16630)
-- Name: pdb_entry_x_emd_entry pdb_entry_x_emd_entry_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pdb_entry_x_emd_entry
    ADD CONSTRAINT pdb_entry_x_emd_entry_pkey PRIMARY KEY (emd_entry_id, pdb_entry_id);


--
-- TOC entry 2796 (class 2606 OID 16575)
-- Name: representation representation_name_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.representation
    ADD CONSTRAINT representation_name_key UNIQUE (name);


--
-- TOC entry 2798 (class 2606 OID 16573)
-- Name: representation representation_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.representation
    ADD CONSTRAINT representation_pkey PRIMARY KEY (id);


--
-- TOC entry 2828 (class 2606 OID 16886)
-- Name: search_history search_history_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.search_history
    ADD CONSTRAINT search_history_pkey PRIMARY KEY (id);


--
-- TOC entry 2808 (class 2606 OID 16648)
-- Name: time_stamp time_stamp_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.time_stamp
    ADD CONSTRAINT time_stamp_pkey PRIMARY KEY (emd_entry_id);


--
-- TOC entry 2788 (class 2606 OID 16531)
-- Name: type_descriptor type_descriptor_name_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.type_descriptor
    ADD CONSTRAINT type_descriptor_name_key UNIQUE (name);


--
-- TOC entry 2790 (class 2606 OID 16529)
-- Name: type_descriptor type_descriptor_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.type_descriptor
    ADD CONSTRAINT type_descriptor_pkey PRIMARY KEY (id);


--
-- TOC entry 2810 (class 2606 OID 16658)
-- Name: update update_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.update
    ADD CONSTRAINT update_pkey PRIMARY KEY (last_update);


--
-- TOC entry 2820 (class 2606 OID 16788)
-- Name: user2 user2_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.user2
    ADD CONSTRAINT user2_pkey PRIMARY KEY (id);


--
-- TOC entry 2822 (class 2606 OID 16818)
-- Name: user user_email_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public."user"
    ADD CONSTRAINT user_email_key UNIQUE (email);


--
-- TOC entry 2824 (class 2606 OID 16816)
-- Name: user user_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public."user"
    ADD CONSTRAINT user_pkey PRIMARY KEY (id);


--
-- TOC entry 2816 (class 2606 OID 16765)
-- Name: user_role user_role_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.user_role
    ADD CONSTRAINT user_role_pkey PRIMARY KEY (id);


--
-- TOC entry 2814 (class 2606 OID 16733)
-- Name: users_roles users_roles_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.users_roles
    ADD CONSTRAINT users_roles_pkey PRIMARY KEY (id);


--
-- TOC entry 2818 (class 2606 OID 16774)
-- Name: userss userss_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.userss
    ADD CONSTRAINT userss_pkey PRIMARY KEY (id);


--
-- TOC entry 2800 (class 2606 OID 16588)
-- Name: volume_filter volume_filter_name_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.volume_filter
    ADD CONSTRAINT volume_filter_name_key UNIQUE (name);


--
-- TOC entry 2802 (class 2606 OID 16586)
-- Name: volume_filter volume_filter_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.volume_filter
    ADD CONSTRAINT volume_filter_pkey PRIMARY KEY (id);


--
-- TOC entry 2835 (class 2606 OID 16713)
-- Name: benchmark_history benchmark_history_representation_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.benchmark_history
    ADD CONSTRAINT benchmark_history_representation_id_fkey FOREIGN KEY (representation_id) REFERENCES public.representation(id);


--
-- TOC entry 2836 (class 2606 OID 16718)
-- Name: benchmark_history benchmark_history_volume_filter_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.benchmark_history
    ADD CONSTRAINT benchmark_history_volume_filter_id_fkey FOREIGN KEY (volume_filter_id) REFERENCES public.volume_filter(id);


--
-- TOC entry 2841 (class 2606 OID 16866)
-- Name: benchmarks_history benchmarks_history_representation_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.benchmarks_history
    ADD CONSTRAINT benchmarks_history_representation_id_fkey FOREIGN KEY (representation_id) REFERENCES public.representation(id);


--
-- TOC entry 2840 (class 2606 OID 16861)
-- Name: benchmarks_history benchmarks_history_user_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.benchmarks_history
    ADD CONSTRAINT benchmarks_history_user_id_fkey FOREIGN KEY (user_id) REFERENCES public."user"(id);


--
-- TOC entry 2842 (class 2606 OID 16871)
-- Name: benchmarks_history benchmarks_history_volume_filter_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.benchmarks_history
    ADD CONSTRAINT benchmarks_history_volume_filter_id_fkey FOREIGN KEY (volume_filter_id) REFERENCES public.volume_filter(id);


--
-- TOC entry 2830 (class 2606 OID 16553)
-- Name: descriptor descriptor_emd_entry_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.descriptor
    ADD CONSTRAINT descriptor_emd_entry_id_fkey FOREIGN KEY (emd_entry_id) REFERENCES public.emd_entry(id);


--
-- TOC entry 2831 (class 2606 OID 16558)
-- Name: descriptor descriptor_type_descriptor_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.descriptor
    ADD CONSTRAINT descriptor_type_descriptor_id_fkey FOREIGN KEY (type_descriptor_id) REFERENCES public.type_descriptor(id);


--
-- TOC entry 2829 (class 2606 OID 16540)
-- Name: emd_entry emd_entry_map_information_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.emd_entry
    ADD CONSTRAINT emd_entry_map_information_id_fkey FOREIGN KEY (map_information_id) REFERENCES public.map_information(id);


--
-- TOC entry 2833 (class 2606 OID 16636)
-- Name: pdb_entry_x_emd_entry pdb_entry_x_emd_entry_emd_entry_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pdb_entry_x_emd_entry
    ADD CONSTRAINT pdb_entry_x_emd_entry_emd_entry_id_fkey FOREIGN KEY (emd_entry_id) REFERENCES public.emd_entry(id);


--
-- TOC entry 2832 (class 2606 OID 16631)
-- Name: pdb_entry_x_emd_entry pdb_entry_x_emd_entry_pdb_entry_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.pdb_entry_x_emd_entry
    ADD CONSTRAINT pdb_entry_x_emd_entry_pdb_entry_id_fkey FOREIGN KEY (pdb_entry_id) REFERENCES public.pdb_entry(id);


--
-- TOC entry 2844 (class 2606 OID 16892)
-- Name: search_history search_history_emd_entry_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.search_history
    ADD CONSTRAINT search_history_emd_entry_id_fkey FOREIGN KEY (emd_entry_id) REFERENCES public.emd_entry(id);


--
-- TOC entry 2845 (class 2606 OID 16897)
-- Name: search_history search_history_representation_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.search_history
    ADD CONSTRAINT search_history_representation_id_fkey FOREIGN KEY (representation_id) REFERENCES public.representation(id);


--
-- TOC entry 2843 (class 2606 OID 16887)
-- Name: search_history search_history_user_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.search_history
    ADD CONSTRAINT search_history_user_id_fkey FOREIGN KEY (user_id) REFERENCES public."user"(id);


--
-- TOC entry 2846 (class 2606 OID 16902)
-- Name: search_history search_history_volume_filter_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.search_history
    ADD CONSTRAINT search_history_volume_filter_id_fkey FOREIGN KEY (volume_filter_id) REFERENCES public.volume_filter(id);


--
-- TOC entry 2834 (class 2606 OID 16649)
-- Name: time_stamp time_stamp_emd_entry_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.time_stamp
    ADD CONSTRAINT time_stamp_emd_entry_id_fkey FOREIGN KEY (emd_entry_id) REFERENCES public.emd_entry(id);


--
-- TOC entry 2838 (class 2606 OID 16789)
-- Name: user2 user2_role_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.user2
    ADD CONSTRAINT user2_role_fkey FOREIGN KEY (role) REFERENCES public.user_role(id);


--
-- TOC entry 2839 (class 2606 OID 16819)
-- Name: user user_role_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public."user"
    ADD CONSTRAINT user_role_fkey FOREIGN KEY (role) REFERENCES public.user_role(id);


--
-- TOC entry 2837 (class 2606 OID 16775)
-- Name: userss userss_role_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.userss
    ADD CONSTRAINT userss_role_fkey FOREIGN KEY (role) REFERENCES public.user_role(id);


-- Completed on 2019-05-21 14:18:52

--
-- PostgreSQL database dump complete
--

