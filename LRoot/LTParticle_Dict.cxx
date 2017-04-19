// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LTParticle_Dict

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
#include "LTParticle.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_LTParticle(void *p = 0);
   static void *newArray_LTParticle(Long_t size, void *p);
   static void delete_LTParticle(void *p);
   static void deleteArray_LTParticle(void *p);
   static void destruct_LTParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LTParticle*)
   {
      ::LTParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LTParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LTParticle", ::LTParticle::Class_Version(), "LTParticle.hxx", 16,
                  typeid(::LTParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LTParticle::Dictionary, isa_proxy, 4,
                  sizeof(::LTParticle) );
      instance.SetNew(&new_LTParticle);
      instance.SetNewArray(&newArray_LTParticle);
      instance.SetDelete(&delete_LTParticle);
      instance.SetDeleteArray(&deleteArray_LTParticle);
      instance.SetDestructor(&destruct_LTParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LTParticle*)
   {
      return GenerateInitInstanceLocal((::LTParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::LTParticle*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr LTParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LTParticle::Class_Name()
{
   return "LTParticle";
}

//______________________________________________________________________________
const char *LTParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LTParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LTParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LTParticle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTParticle*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void LTParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class LTParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LTParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(LTParticle::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LTParticle(void *p) {
      return  p ? new(p) ::LTParticle : new ::LTParticle;
   }
   static void *newArray_LTParticle(Long_t nElements, void *p) {
      return p ? new(p) ::LTParticle[nElements] : new ::LTParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_LTParticle(void *p) {
      delete ((::LTParticle*)p);
   }
   static void deleteArray_LTParticle(void *p) {
      delete [] ((::LTParticle*)p);
   }
   static void destruct_LTParticle(void *p) {
      typedef ::LTParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LTParticle

namespace {
  void TriggerDictionaryInitialization_LTParticle_Dict_Impl() {
    static const char* headers[] = {
"LTParticle.hxx",
0
    };
    static const char* includePaths[] = {
"/usr/local/root/include",
"/.",
"/usr/local/root-6.08.00-build/include",
"/home/gold/doubleBeta/LRoot/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LTParticle_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$LTParticle.hxx")))  LTParticle;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LTParticle_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "LTParticle.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"LTParticle", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LTParticle_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LTParticle_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LTParticle_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LTParticle_Dict() {
  TriggerDictionaryInitialization_LTParticle_Dict_Impl();
}
