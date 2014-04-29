using System.Web.Mvc;

namespace PhdThesisIngaSamonenko.Controllers
{
    public class BDController : Controller
    {
        public ActionResult Index()
        {
            return File("BD.R", "text/plain");
        }

    }
}
