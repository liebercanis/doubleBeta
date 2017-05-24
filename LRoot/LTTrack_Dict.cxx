// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LTTrack_Dict

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
#include "LTTrack.hxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_LTTrack(void *p = 0);
   static void *newArray_LTTrack(Long_t size, void *p);
   static void delete_LTTrack(void *p);
   static void deleteArray_LTTrack(void *p);
   static void destruct_LTTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LTTrack*)
   {
      ::LTTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LTTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LTTrack", ::LTTrack::Class_Version(), "LTTrack.hxx", 22,
                  typeid(::LTTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LTTrack::Dictionary, isa_proxy, 4,
                  sizeof(::LTTrack) );
      instance.SetNew(&new_LTTrack);
      instance.SetNewArray(&newArray_LTTrack);
      instance.SetDelete(&delete_LTTrack);
      instance.SetDeleteArray(&deleteArray_LTTrack);
      instance.SetDestructor(&destruct_LTTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LTTrack*)
   {
      return GenerateInitInstanceLocal((::LTTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::LTTrack*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr LTTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LTTrack::Class_Name()
{
   return "LTTrack";
}

//______________________________________________________________________________
const char *LTTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LTTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LTTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LTTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LTTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LTTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void LTTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class LTTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LTTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(LTTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LTTrack(void *p) {
      return  p ? new(p) ::LTTrack : new ::LTTrack;
   }
   static void *newArray_LTTrack(Long_t nElements, void *p) {
      return p ? new(p) ::LTTrack[nElements] : new ::LTTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_LTTrack(void *p) {
      delete ((::LTTrack*)p);
   }
   static void deleteArray_LTTrack(void *p) {
      delete [] ((::LTTrack*)p);
   }
   static void destruct_LTTrack(void *p) {
      typedef ::LTTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LTTrack

namespace {
  void TriggerDictionaryInitialization_LTTrack_Dict_Impl() {
    static const char* headers[] = {
"LTTrack.hxx",
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
#line 1 "LTTrack_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$LTTrack.hxx")))  LTTrack;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LTTrack_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "LTTrack.hxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"LTTrack", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LTTrack_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LTTrack_Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LTTrack_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LTTrack_Dict() {
  TriggerDictionaryInitialization_LTTrack_Dict_Impl();
}
