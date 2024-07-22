import itertools
import logging
import os
import re
import string
import sys
import urllib
import urllib.error
import urllib.request
from typing import Optional
from time import sleep

import requests
from bs4 import BeautifulSoup
import bs4

from cazy_parser.utils import init_bar

log = logging.getLogger("cazylog")


def process_ec_text(ec_text: str) -> str:
    ec_re = re.compile(r"(\d{1}.\d+.[0-9a-zA-Z-]+.[0-9a-zA-Z-]+(?!\.))")
    ec_matches = ec_re.findall(ec_text)

    return ":".join(ec_matches)

def process_pdb_text(pdb_text: str) -> str:
    pdb_red = re.compile(r"(\w{4,8})\[(.+?)\]")

    pdb_chains = [
        f"{pdb_id}_{chain}"
        for pdb_id, chains in pdb_red.findall(pdb_text)
            for chain in chains.split(",")
    ]

    return ":".join(pdb_chains)

def process_uniprot_text(uniprot_text: str) -> str:
    """TODO make sure this works across all cases"""
    processed_ids = []
    uniprot_text = uniprot_text.replace("\xa0","")
    while uniprot_text != "":
        if uniprot_text.startswith("A0"):
            processed_ids.append(uniprot_text[:10])
            uniprot_text = uniprot_text[10:]
        elif len(uniprot_text) % 6 == 0:
            processed_ids.append(uniprot_text[:6])
            uniprot_text = uniprot_text[6:]
        else:
            log.warning(f"Uniprot text {uniprot_text} is not divisble by 6; not sure if correct")
            processed_ids.append(uniprot_text)
            uniprot_text = ""


    return ":".join(processed_ids)

def fetch_data(link_list):
    """
    Parse list of links and extract information.

    Parameters
    ----------
    link_list : list
        List of links to parse


    Returns
    -------
    data : list
        List of dictionaries containing information about each entry.

    """

    data = []
    for _, link in enumerate(link_list):
        if ".txt" in link:
            link_data = get_data_from_txt(link)
        else:
            link_data = get_data_from_html(link)

        data.append(link_data)

    return list(itertools.chain(*data))


def parse_td(td_list):
    """
    Parse a table from the HTML page.

    Parameters
    ----------
    td_list : list
        List of table elements.

    Returns
    -------
    data_dict : dic
        Dictionary of data.

    """
    data_dict = {}
    if len(td_list) <= 1 or td_list[0].text == "Protein Name":
        # this is likely the header
        return data_dict

    protein_name = td_list[0].text
    ec = process_ec_text(td_list[1].text)
    # referece = td_list[2].text
    organism = td_list[3].text
    try:
        genbank = td_list[4].find("a").text
    except AttributeError:
        genbank = "unavailable"

    uniprot = process_uniprot_text(td_list[5].text)
    pdb = process_pdb_text(td_list[6].text)

    subfamily = td_list[7].text if len(td_list) > 7 else "unavailable"

    data_dict["protein_name"] = protein_name
    data_dict["ec"] = ec
    data_dict["organism"] = organism
    data_dict["genbank"] = genbank
    data_dict["uniprot"] = uniprot
    data_dict["pdb"] = pdb
    data_dict["subfamily"] = subfamily

    return data_dict


def parse_table(table: bs4.element.Tag) -> list[dict]:
    """
    Parse a beautiful soup table and retrieve information.

    Parameters
    ----------
    table : bs4.element.Tag
        Beautiful soup table.

    Returns
    -------
    table_data : list
        List of dictionaries containing information from the table.

    """
    table_data = []
    for tr in table.findAll("tr"):
        tds = tr.findAll("td")
        td_dic = parse_td(tds)
        table_data.append(td_dic)
    return table_data


def get_data_from_html(link):
    """
    Retrieve information from the HTML page.

    Parameters
    ----------
    link : str
        Link to the page.

    """
    soup = BeautifulSoup(urllib.request.urlopen(link), features="html.parser")
    table = soup.find("table", attrs={"class": "listing"})
    domain = ""
    family = link.split(".org/")[-1].split("_")[0]
    data_list = parse_table(table)

    tag = "characterized" if "characterized" in link else ""
    # add more information to the data
    for data in data_list:
        data["tag"] = tag
        data["family"] = family
        data["domain"] = domain

    return data_list


def get_data_from_txt(link):
    """
    Retrieve information from the TXT file.

    Parameters
    ----------
    link : str
        Link to the page.

    Returns
    -------
    data_list : list
        List of dictionaries containing information from the TXT file.

    """
    data_list = []
    response = requests.get(link)
    tag = "characterized" if "characterized" in link else ""
    for line in response.text.split(os.linesep):
        data_dict = {}
        data = line.split("\t")

        try:
            family = data[0]
            domain = data[1]
            species = data[2]
            gene = data[3]

            data_dict["family"] = family
            data_dict["domain"] = domain
            data_dict["species"] = species
            data_dict["genbank"] = gene
            data_dict["tag"] = tag
        except IndexError:
            # no data for this line
            continue

        data_list.append(data_dict)

    return data_list


def fetch_links(
    enzyme_class: str,
    family: Optional[int] = None,
    subfamily: Optional[int] = None,
    characterized: bool = False,
    **kwargs
) -> list[str]:
    """Fetch link structure for an enzyme class."""

    main_class_link = f"http://www.cazy.org/{enzyme_class}.html"
    log.info(f"Fetching links for {enzyme_class}, url: {main_class_link}")
    family_list = fetch_families(main_class_link)


    # Filter by family
    if family:
        if subfamily:
            log.info(f"Only using links of family {family} subfamily {subfamily}")
            family_list = [e for e in family_list if e[2:] == f"{family}_{subfamily}"]
        else:
            log.info(f"Only using links of family {family}")
            family_list = [e for e in family_list if int(e[2:]) == family]

    log.info(f"Found {len(family_list)} families")

    if characterized:
        log.info("Only using characterized links")

    if not family_list:
        log.error("No links were found.")
        sys.exit()

    sleep_time = kwargs.get("sleep_time", None)
    page_list = []
    try: 
        for j, family in enumerate(family_list, start=1):
            sleep(1)
            if isinstance(sleep_time, int) and j % 10 == 0:
                log.debug(f"Completed {j}: Sleeping for {sleep_time} seconds...")
                sleep(sleep_time)
            # bar.update(j)
            family_link = f"http://www.cazy.org/{family}.html"

            # TODO: Implement checkpoint for link fetching

            family_soup = BeautifulSoup(
                urllib.request.urlopen(family_link), features="html.parser"
            )

            # Find the links to the individual pages
            superfamily_links = []
            for line in family_soup.findAll("span", attrs={"class": "choix"}):
                _link = line.find("a")["href"]
                if "krona" not in _link and "structure" not in _link:
                    superfamily_links.append(_link)

            log.debug(f"Superfamily links: {superfamily_links}")
            for link in superfamily_links:
                page_zero = link
                try:
                    soup = BeautifulSoup(
                        urllib.request.urlopen(link), features="html.parser"
                    )
                except ValueError:
                    # This is a link to a text file
                    page_list.append(f"http://www.cazy.org/{link}")
                    continue

                # Get page list for the family
                page_index_list = soup.findAll(name="a", attrs={"class": "lien_pagination"})
                if bool(page_index_list):
                    # =====================#
                    # be careful with this #
                    first_page_idx = int(re.findall(r"=(\d*)#", str(page_index_list[0]))[0])
                    last_page_idx = int(re.findall(r"=(\d*)#", str(page_index_list[-2]))[0])
                    # =====================#

                    page_list.append(page_zero)
                    page_range = range(
                        first_page_idx, last_page_idx + first_page_idx, first_page_idx
                    )
                    for i in page_range:
                        sub_str = page_index_list[0]["href"].split("=")[0]
                        link = f"http://www.cazy.org/{sub_str}={i}"
                        log.debug(f"Page index links: {link}")
                        page_list.append(link)

                else:
                    page_list.append(page_zero)

    except urllib.error.HTTPError as e:
        log.error(f"HTTP Error: {e}")
        log.warning(f"Completed families: {family_list[:j]}")
        log.warning(f"Incomplete families: {family_list[j:]}")
        log.warning(f"Moving forward with completed links: {len(page_list)}")

    if characterized:
        page_list = [e for e in page_list if "characterized" in e]

    log.info(f"Total links found: {len(page_list)}")

    return page_list


def fetch_families(link):
    """
    Identify family link structure and return a list.

    Parameters
    ----------
    link : str
        Link to the page.

    Returns
    -------
    family_link_list : list
        List of family links.

    """
    enzyme_regex = r"([A-Z]{2}\d*_?\d?).html"
    soup = BeautifulSoup(urllib.request.urlopen(link), features="html.parser")
    all_links = soup.find_all("a")
    family_link_list = []
    for link in all_links:
        try:
            href = re.findall(enzyme_regex, link.attrs["href"])[0]
            family_link_list.append(href)
        except (IndexError, KeyError):
            continue

    return list(set(family_link_list))


def fetch_species():
    """Gather species names and IDs (full genome only)."""
    log.info("> Gathering species with full genomes")
    # a = archea // b = bacteria // e = eukaryota // v = virus
    species_domain_list = ["a", "b", "e", "v"]
    species_dic = {}

    bar = init_bar()
    bar.max_value = len(string.ascii_uppercase) * len(species_domain_list)
    bar.start()

    counter = 0
    for domain in species_domain_list:
        for initial in string.ascii_uppercase:
            counter += 1
            bar.update(counter)
            link = f"http://www.cazy.org/{domain}{initial}.html"
            f = urllib.request.urlopen(link)
            species_list_hp = f.read().decode("utf-8")
            # parse webpage
            index_list = re.findall(
                rf'"http://www.cazy.org/{species_domain_list}(\d.*).html"'
                r' class="nav">(.*)</a>',
                species_list_hp,
            )
            for entry in index_list:
                index, species = entry
                try:
                    species_dic[species].append(index)
                except KeyError:
                    species_dic[species] = [index]
    bar.finish()

    # Double check to see which of the species codes are valid
    for j, species in enumerate(species_dic.keys()):
        entry_list = species_dic[species]
        if len(entry_list) > 1:
            # More than one entry for this species
            #  > This is (likely) a duplicate
            #  > Use the higher number, should be the newer page
            newer_entry = max([int(i.split("b")[-1]) for i in entry_list])
            selected_entry = "b%i" % newer_entry

            species_dic[species] = selected_entry
        else:
            species_dic[species] = species_dic[species][0]

    return species_dic


def retrieve_cazy_metadata(
    enzyme_name: str, 
    family: Optional[int], 
    subfamily: Optional[int], 
    characterized: bool = False,
    filter_: Optional[str] = None,
    **kwargs
) -> list[str]:
    """Retrieve metadata of given enzyme(s)."""
    page_list = fetch_links(enzyme_name, family, subfamily, characterized, **kwargs)
    data = fetch_data(page_list)
    log.info(f"Retrieved {len(data)} entries")

    if filter_:
        data = [element for element in data if filter_ in element]
        log.info(f"Retrieved {len(data)} entries with valid {filter_}")

    return data

def retrieve_prop_from_metadata(
        property: str, 
        enzyme_name: str,  
        family: Optional[int], 
        subfamily: Optional[int], 
        characterized: bool = False
) -> list:
    """Retrieve property for a given enzyme(s)."""
    metadata = retrieve_cazy_metadata(enzyme_name, family, subfamily, characterized, property)
    
    prop_data = [
        element[property]
        for element in metadata
        if property in element
    ]
    
    log.info(f"Retrieved {len(prop_data)} entries with valid {property}")

    return prop_data

def retrieve_uniprot_from_metadata(
    enzyme_name: str, family: Optional[int], subfamily: Optional[int], characterized: bool = False
) -> list[str]:
    """Retrieve uniprot IDs for a given enzyme(s)."""
    return retrieve_prop_from_metadata("uniprot", enzyme_name, family, subfamily, characterized)
    
def retrieve_genbank_from_metadata(
    enzyme_name: str, family: Optional[int], subfamily: Optional[int], characterized: bool = False
) -> list[str]:
    """Retrieve uniprot IDs for a given enzyme(s)."""
    return retrieve_prop_from_metadata("genbank", enzyme_name, family, subfamily, characterized)
    

# def retrieve_metadata_if_genbank_id(
#     enzyme_name: str, family: Optional[int], subfamily: Optional[int], characterized: bool = False
# ) -> list[str]:
#     """Retrieve genbank IDs for a given enzyme."""
#     page_list = fetch_links(enzyme_name, family, subfamily, characterized)
#     data = fetch_data(page_list)
#     data_list = [element for element in data if "genbank" in element]

#     log.info(f"Retrieved {len(data_list)} entries with valid genbank ids")

#     return data_list



    
