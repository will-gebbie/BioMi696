#!/usr/bin/env python3

from serpapi import GoogleSearch

name1 = 'Kalyuzhnaya'
name2 = 'Collins'

params = {
    "engine": "google_scholar",
    "q": "author:{} author:{}".format(name1, name2),
    "as_vis": 1,
    "hl": "en",
    "api_key": "Secret!!"
}

f = open('{}_{}.txt'.format(name1,name2), "a")

search = GoogleSearch(params)
results = search.get_dict()

organic = results["organic_results"]


q_citation_list = []
for key in organic:
    links = key['inline_links']
    cite_link = links['serpapi_cite_link']
    q_idx = cite_link.index('&q=')
    q_citation_list.append(cite_link[q_idx+3:])


for q in q_citation_list:
    cite_params = {
        "engine": "google_scholar_cite",
        "q": "{}".format(q),
        "as_vis": 1,
        "hl": "en",
        "api_key": "Secret!"
    }

    cite_search = GoogleSearch(cite_params)
    cite_results = cite_search.get_dict()
    citations = cite_results["citations"]

    print(citations)
    f.write(citations[1]['snippet'])
    f.write('\n\n')


f.close()