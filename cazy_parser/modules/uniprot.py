from io import StringIO
import logging
import re

from Bio import SeqIO
import requests

log = logging.getLogger("cazylog")


UNIPROT_FASTA_URL = "https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
UNISAVE_FASTA_URL = "https://rest.uniprot.org/unisave/{uniprot_id}?format=fasta"
UNIPROT_BATCH_URL= "https://rest.uniprot.org//uniprotkb/accessions?format=fasta&accessions={uniprot_ids}"

def fetch_uniprot_sequence(uniprot_id: str, endpoint=UNIPROT_FASTA_URL) -> tuple[str, str]:
    url = endpoint.format(uniprot_id=uniprot_id)
    response = requests.get(url)
    log.debug(f"{url}: {response.status_code} text: {bool(response.text)}")
    if response.status_code == 200:
        text = response.text
        if not response.text:
            return fetch_uniprot_sequence(uniprot_id, endpoint=UNISAVE_FASTA_URL)
        
        if text.find("\n>",1) != -1:
            text = text[:text.find("\n>",1)]

        # print(text)
        header, sequence = text.split("\n", maxsplit=1)
        return header, sequence.replace("\n", "").strip()
    
    elif response.status_code == 301:
       return fetch_uniprot_sequence(uniprot_id, endpoint=UNISAVE_FASTA_URL)
    else:
        return None, None
    

def batch_fetch_uniprot_sequence(uniprot_ids: list[str]):
    """This fails for IDs that require pinging the UNISAVE endpoint"""
    ids = ",".join(set(uniprot_ids))

    url = UNIPROT_BATCH_URL.format(uniprot_ids=ids)
    response = requests.get(url)

    if response.status_code == 200:
        fasta_batch = StringIO(response.text)

        for record in SeqIO.parse(fasta_batch, "fasta"):
            print(f"ID: {record.id}")
            print(f"Description: {record.description}")
            print(f"Sequence: {record.seq}")


    

# TODO asyncio implementation from chat GPT; still need to test and implement

# import asyncio
# from aiohttp import ClientSession, ClientResponseError, ClientTimeout

# async def fetch_uniprot_sequence(session: ClientSession, uniprot_id: str):
#     url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
#     try:
#         async with session.get(url) as response:
#             response.raise_for_status()
#             text = await response.text()
#             header, sequence = text.split("\n", maxsplit=1)
#             return header, sequence.replace("\n", "").strip()
#     except ClientResponseError as e:
#         print(f"Failed to fetch {uniprot_id}: {e}")
#         return None, None

# async def fetch_all_sequences(uniprot_ids):
#     async with ClientSession(timeout=ClientTimeout(total=60)) as session:
#         tasks = []
#         for uniprot_id in uniprot_ids:
#             tasks.append(fetch_uniprot_sequence(session, uniprot_id))
        
#         # Limit the number of concurrent tasks
#         results = []
#         semaphore = asyncio.Semaphore(10)
        
#         async def sem_task(task):
#             async with semaphore:
#                 return await task
        
#         results = await asyncio.gather(*[sem_task(task) for task in tasks])
#         return results

# # List of UniProt IDs
# uniprot_ids = [
#     'P12345', 'Q8N158', 'P67890',  # Add your IDs here
# ]

# # Fetch and save sequences asynchronously
# async def main():
#     results = await fetch_all_sequences(uniprot_ids)
#     with open('sequences.fasta', 'w') as out_file:
#         for header, sequence in results:
#             if header and sequence:
#                 out_file.write(f"{header}\n{sequence}\n")

# if __name__ == "__main__":
#     asyncio.run(main())
