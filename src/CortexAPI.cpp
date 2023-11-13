#include "../inc/CortexAPI.h"

size_t CortexAPI::WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

CortexAPI::CortexAPI() {
    curl_global_init(CURL_GLOBAL_DEFAULT);
}

CortexAPI::~CortexAPI() {
    curl_global_cleanup();
}

std::vector<NType> CortexAPI::postImage(const std::string& imagePath) {
    CURL* curl;
    CURLcode res;
    std::string readBuffer;

    struct curl_httppost* post = NULL;
    struct curl_httppost* last = NULL;

    curl = curl_easy_init();

    if(curl) {
        curl_formadd(&post, &last,
                     CURLFORM_COPYNAME, "image",
                     CURLFORM_FILE, imagePath.c_str(),
                     CURLFORM_CONTENTTYPE, "image/jpeg",
                     CURLFORM_END);

        curl_easy_setopt(curl, CURLOPT_URL, "https://7m15gatms9.execute-api.us-east-1.amazonaws.com/v1/embeddings");
        curl_easy_setopt(curl, CURLOPT_HTTPPOST, post);
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L);
        struct curl_slist* list = NULL;
        list = curl_slist_append(list, "x-api-key: eJLlLzKise4B1Bl7DC6XQrbGcIpShDD5eaC9XNl2");
        curl_easy_setopt(curl, CURLOPT_HTTPHEADER, list);

        res = curl_easy_perform(curl);
        if(res != CURLE_OK) {
            std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << std::endl;
        }

        curl_easy_cleanup(curl);
        curl_formfree(post);
        curl_slist_free_all(list);
    }

    std::stringstream ss(readBuffer);
    std::string result = ss.str();
    std::istringstream iss(result);
    iss.get();
    std::string token;
    std::vector<NType> float_values;

    while (std::getline(iss, token, ',')) {
        if (!token.empty()) {
            float value;
            if (std::istringstream(token) >> value) {
                float_values.emplace_back(value);
            }
        }
    }

    return float_values;
}