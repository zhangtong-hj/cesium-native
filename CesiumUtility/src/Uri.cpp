#include "CesiumUtility/Uri.h"

#include <uriparser/Uri.h>

#include <stdexcept>

bool isChinese(char c) {
  int x = c & 0xE0;
  if (x == 224)
    return true;
  else
    return false;
}

 static char _HEX_[] = {
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0x2b, 0x21, 0,    0x23,
    0x24, 0,    0x26, 0x27, 0x28, 0x29, 0x2a, 0,    0x2c, 0x2d, 0x2e, 0x2f,
    0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3a, 0x3b,
    0,    0x3d, 0,    0x3f, 0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47,
    0x48, 0x49, 0x4a, 0x4b, 0x4c, 0x4d, 0x4e, 0x4f, 0x50, 0x51, 0x52, 0x53,
    0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5a, 0x5b, 0,    0x5d, 0,    0x5f,
    0,    0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x6a, 0x6b,
    0x6c, 0x6d, 0x6e, 0x6f, 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77,
    0x78, 0x79, 0x7a, 0,    0,    0,    0x7e, 0}; // RFC3986
std::string EncodeURL(const std::string& s) {
   std::string r;
   int location = 0;
   bool isFind = false;
   for (char c : s) {
     if (c == '%' && (location + 2) < s.size()) {
       if (s[location + 1] == '2' && s[location + 2] == 'B') {
         isFind = true;
       }
     }
     if (c == '[') {
       r.push_back('%');
       r.push_back('5');
       r.push_back('b');
       isFind = true;
     }
     if (c == ']') {
       r.push_back('%');
       r.push_back('5');
       r.push_back('d');
       isFind = true;
     }
     if (!isFind) {
       if (c > '\377') {
         if (_HEX_[c])
           r.push_back(_HEX_[c]); // if (_RFC2396[c])r.push_back(c);
         else {
           r.push_back(0x25);
           char o = (c & 0xF0) >> 4;
           o += o > 9 ? 0x37 : 0x30;
           r.push_back(o);
           o = c & 0x0F;
           o += o > 9 ? 0x37 : 0x30;
           r.push_back(o);
         }
       } else {
         r.push_back(0x25);
         char o = (static_cast<uint8_t>(c) & 0xF0) >> 4;
         o += o > 9 ? 0x37 : 0x30;
         r.push_back(o);
         o = static_cast<uint8_t>(c) & 0x0F;
         o += o > 9 ? 0x37 : 0x30;
         r.push_back(o);
       }
     } else {
       if (c != '[' && c !=']') {
         r.push_back(c);
       }
     }

     location++;
   }
   return r;
 }


namespace CesiumUtility {
std::string Uri::resolve(
    const std::string& base,
    const std::string& relative,
    bool useBaseQuery) {
  UriUriA baseUri;

  const char* errorPos = NULL;
  auto EncodedBaseUri = EncodeURL(base);
  if (uriParseSingleUriA(&baseUri, EncodedBaseUri.c_str(), &errorPos) !=
      URI_SUCCESS) {
    // Could not parse the base, so just use the relative directly and hope for
    // the best.
    return relative;
  }

  auto Encodedrelative = EncodeURL(relative);
  UriUriA relativeUri;
  if (uriParseSingleUriA(&relativeUri, Encodedrelative.c_str(), nullptr) !=
      URI_SUCCESS) {
    // Could not parse one of the URLs, so just use the relative directly and
    // hope for the best.
    uriFreeUriMembersA(&baseUri);
    return relative;
  }

  UriUriA resolvedUri;
  if (uriAddBaseUriA(&resolvedUri, &relativeUri, &baseUri) != URI_SUCCESS) {
    uriFreeUriMembersA(&resolvedUri);
    uriFreeUriMembersA(&relativeUri);
    uriFreeUriMembersA(&baseUri);
    return relative;
  }

  if (uriNormalizeSyntaxA(&resolvedUri) != URI_SUCCESS) {
    uriFreeUriMembersA(&resolvedUri);
    uriFreeUriMembersA(&relativeUri);
    uriFreeUriMembersA(&baseUri);
    return relative;
  }

  int charsRequired;
  if (uriToStringCharsRequiredA(&resolvedUri, &charsRequired) != URI_SUCCESS) {
    uriFreeUriMembersA(&resolvedUri);
    uriFreeUriMembersA(&relativeUri);
    uriFreeUriMembersA(&baseUri);
    return relative;
  }

  std::string result(static_cast<size_t>(charsRequired), ' ');

  if (uriToStringA(
          const_cast<char*>(result.c_str()),
          &resolvedUri,
          charsRequired + 1,
          nullptr) != URI_SUCCESS) {
    uriFreeUriMembersA(&resolvedUri);
    uriFreeUriMembersA(&relativeUri);
    uriFreeUriMembersA(&baseUri);
    return relative;
  }

  if (useBaseQuery) {
    std::string query(baseUri.query.first, baseUri.query.afterLast);
    if (query.length() > 0) {
      if (resolvedUri.query.first) {
        result += "&" + query;
      } else {
        result += "?" + query;
      }
    }
  }

  uriFreeUriMembersA(&resolvedUri);
  uriFreeUriMembersA(&relativeUri);
  uriFreeUriMembersA(&baseUri);

  return result;
}

std::string Uri::addQuery(
    const std::string& uri,
    const std::string& key,
    const std::string& value) {
  // TODO
  if (uri.find('?') != std::string::npos) {
    return uri + "&" + key + "=" + value;
  }
  return uri + "?" + key + "=" + value;
  // UriUriA baseUri;

  // if (uriParseSingleUriA(&baseUri, uri.c_str(), nullptr) != URI_SUCCESS)
  //{
  //	// TODO: report error
  //	return uri;
  //}

  // uriFreeUriMembersA(&baseUri);
}

std::string Uri::substituteTemplateParameters(
    const std::string& templateUri,
    const std::function<SubstitutionCallbackSignature>& substitutionCallback) {
  std::string result;
  std::string placeholder;

  size_t startPos = 0;
  size_t nextPos;

  // Find the start of a parameter
  while ((nextPos = templateUri.find('{', startPos)) != std::string::npos) {
    result.append(templateUri, startPos, nextPos - startPos);

    // Find the end of this parameter
    ++nextPos;
    const size_t endPos = templateUri.find('}', nextPos);
    if (endPos == std::string::npos) {
      throw std::runtime_error("Unclosed template parameter");
    }

    placeholder = templateUri.substr(nextPos, endPos - nextPos);
    result.append(substitutionCallback(placeholder));

    startPos = endPos + 1;
  }

  result.append(templateUri, startPos, templateUri.length() - startPos);

  return result;
}

std::string Uri::escape(const std::string& s) {
  // In the worst case, escaping causes each character to turn into three.
  std::string result(s.size() * 3, '\0');
  char* pTerminator = uriEscapeExA(
      s.data(),
      s.data() + s.size(),
      result.data(),
      URI_FALSE,
      URI_FALSE);
  result.resize(size_t(pTerminator - result.data()));
  return result;
}
} // namespace CesiumUtility
