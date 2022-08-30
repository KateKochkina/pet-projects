from elasticsearch import Elasticsearch

es = Elasticsearch()
doc_content = set()
count_repeating = 0

res = es.search(size=10000)
print('total: ', res['hits']['total'])
for hit in res['hits']['hits']:
    es.delete(index=hit['_index'], doc_type=hit['_type'], id=hit['_id'])
    # content = hit['_source']['page_content']
    # if content in doc_content:
    #     count_repeating += 1
    #     print(hit['_source']['url'])
    # doc_content.add(content)

#print('count_repeating: ', count_repeating)
