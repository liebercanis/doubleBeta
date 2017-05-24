// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LTPVertex_Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "LTPVertex.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_LTPVertex(void *p = 0);
   static void *newArray_LTPVertex(Long_t size, void *p);
   static void delete_LTPVertex(void *p);
   static void deleteArray_LTPVertex(void *p);
   static void destruct_LTPVertex(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LTPVertex*)
   {
      ::LTPVertex *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LTPVertex >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LTPVertex", ::LTPVertex::Class_Version(), "LTPVertex.hxx", 12,
                  typeid(::LTPVertex), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LTPVertex::Dictionary, isa_proxy, 4,
                  sizeof(::LTPVertex) );
      instance.SetNew(&new_LTPVertex);
      instance.SetNewArray(&newArray_LTPVertex);
      instance.SetDelete(&delete_LTPVertex);
      instance.SetDeleteArray(&deleteArray_LTPVertex);
      instance.SetDestructor(&destruct_LTPVertex);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LTPVertex*)
   {
      return GenerateInitInstanceLocal((::LTPVertex*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::LTPVertex*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr LTPVertex::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LTPVertex::Class_Name()
{
   return "LTPVertex";
}

//______________________________________________________________________________
const char *LTPVertex::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTPVertex*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LTPVertex::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTPVertex*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LTPVertex::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTPVertex*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LTPVertex::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTPVertex*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void LTPVertex::Streamer(TBuffer &R__b)
{
   // Stream an object of class LTPVertex.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LTPVertex::Class(),this);
   } else {
      R__b.WriteClassBuffer(LTPVertex::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LTPVertex(void *p) {
      return  p ? new(p) ::LTPVertex : new ::LTPVertex;
   }
   static void *newArray_LTPVertex(Long_t nElements, void *p) {
      return p ? new(p) ::LTPVertex[nElements] : new ::LTPVertex[nElements];
   }
   // Wrapper around operator delete
   static void delete_LTPVertex(void *p) {
      delete ((::LTPVertex*)p);
   }
   static void deleteArray_LTPVertex(void *p) {
      delete [] ((::LTPVertex*)p);
   }
   static void destruct_LTPVertex(void *p) {
      typedef ::LTPVertex current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LTPVertex

namespace ROOT {
   static TClass *vectorlELTParticlegR_Dictionary();
   static void vectorlELTParticlegR_TClassManip(TClass*);
   static void *new_vectorlELTParticlegR(void *p = 0);
   static void *newArray_vectorlELTParticlegR(Long_t size, void *p);
   static void delete_vectorlELTParticlegR(void *p);
   static void deleteArray_vectorlELTParticlegR(void *p);
   static void destruct_vectorlELTParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LTParticle>*)
   {
      vector<LTParticle> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LTParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LTParticle>", -2, "vector", 210,
                  typeid(vector<LTParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELTParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<LTParticle>) );
      instance.SetNew(&new_vectorlELTParticlegR);
      instance.SetNewArray(&newArray_vectorlELTParticlegR);
      instance.SetDelete(&delete_vectorlELTParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlELTParticlegR);
      instance.SetDestructor(&destruct_vectorlELTParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LTParticle> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<LTParticle>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELTParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<LTParticle>*)0x0)->GetClass();
      vectorlELTParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELTParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELTParticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LTParticle> : new vector<LTParticle>;
   }
   static void *newArray_vectorlELTParticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<LTParticle>[nElements] : new vector<LTParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELTParticlegR(void *p) {
      delete ((vector<LTParticle>*)p);
   }
   static void deleteArray_vectorlELTParticlegR(void *p) {
      delete [] ((vector<LTParticle>*)p);
   }
   static void destruct_vectorlELTParticlegR(void *p) {
      typedef vector<LTParticle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<LTParticle>

namespace {
  void TriggerDictionaryInitialization_LTPVertex_Dict_Impl() {
    static const char* headers[] = {
"LTPVertex.hxx",
0
    };
    static const char* includePaths[] = {
"/usr/local/root/include",
"/.",
"/usr/local/root-6.08.00-build/include",
"/home/nmcfadde/doubleBeta/LRoot/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LTPVertex_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$LTPVertex.hxx")))  LTParticle;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$LTPVertex.hxx")))  LTPVertex;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LTPVertex_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "LTPVertex.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"LTPVertex", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LTPVertex_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LTPVertex_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LTPVertex_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LTPVertex_Dict() {
  TriggerDictionaryInitialization_LTPVertex_Dict_Impl();
}
